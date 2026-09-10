#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <memory>
#include <raylib.h>
#include <stdint.h>
#include <unistd.h>

constexpr std::size_t SIZE = 300;
constexpr std::size_t N = SIZE * SIZE;
constexpr float cs2 = 1.0f / 3.0f;

struct Vector {
  float x, y;
  Vector(float _x, float _y) : x{_x}, y{_y} {}
  Vector() : x{0}, y{0} {}
};

// D2Q9 lattice constants as flat arrays indexed by direction 0..8
namespace D2Q9 {
static constexpr std::array<float, 9> ex = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
static constexpr std::array<float, 9> ey = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
static constexpr std::array<float, 9> w = {1.0f / 36, 1.0f / 9, 1.0f / 36,
                                           1.0f / 9,  4.0f / 9, 1.0f / 9,
                                           1.0f / 36, 1.0f / 9, 1.0f / 36};

static constexpr int dir_index(int dx, int dy) {
  return (dx + 1) * 3 + (dy + 1);
}
} // namespace D2Q9

class Lattice {
public:
  // Changes AoS to SoA to make it SIMD vectorizable.
  std::array<std::array<float, N>, 9> old_f;
  std::array<std::array<float, N>, 9> new_f;

  std::array<float, N> density;
  std::array<float, N> vel_x;
  std::array<float, N> vel_y;

  float relaxation_time;
  float min_density, max_density, min_speed, max_speed;

  static std::size_t idx(int x, int y) {
    return static_cast<std::size_t>(x) * SIZE + static_cast<std::size_t>(y);
  }

  // The peak:ambient ratio is a pressure ratio, and BGK is only valid while
  // the flow it drives stays well below the lattice sound speed (cs =
  // 0.577).
  Lattice(float _relaxation_time, float _peak_density = 0.2f,
          float _base_density = 0.15f, float _blob_radius = SIZE * 0.15f,
          Vector _mean_velocity = Vector(0.0f, 0.0f))
      : relaxation_time{_relaxation_time}, min_density{FLT_MAX},
        max_density(-FLT_MAX), min_speed{FLT_MAX}, max_speed(0.0f) {
    init(_peak_density, _base_density, _blob_radius, _mean_velocity);
  }

  static float compute_equilibrium(int d, float density, float ux, float uy) {
    const float edotu = D2Q9::ex[d] * ux + D2Q9::ey[d] * uy;
    const float u2 = ux * ux + uy * uy;
    return (1 + edotu / cs2 + (edotu * edotu) / (2 * cs2 * cs2) -
            u2 / (2 * cs2)) *
           density * D2Q9::w[d];
  }

  // Starts the fluid at rest as a radially symmetric blob: peak_density at
  // the centre of the domain falling off to ambient_density at the edges.
  // The resulting pressure gradient is what drives the flow on the first
  // steps.
  void init(float peak_density, float base_density, float blob_radius,
            Vector mean_velocity) {
    const float cx = (SIZE - 1) * 0.5f;
    const float cy = (SIZE - 1) * 0.5f;
    const float two_r2 = 2.0f * blob_radius * blob_radius;

    for (int x = 0; x < SIZE; x++) {
      for (int y = 0; y < SIZE; y++) {
        const float dx = x - cx;
        const float dy = y - cy;
        // A Gaussian rather than a hard-edged disc: the discontinuity of a
        // sharp edge launches a ringing pressure wave on the first step (That
        // took a while for me to fix :( )
        const float falloff = std::exp(-(dx * dx + dy * dy) / two_r2);
        const float _density =
            base_density + (peak_density - base_density) * falloff;
        const std::size_t i = idx(x, y);

        for (int d = 0; d < 9; d++) {
          old_f[d][i] = compute_equilibrium(d, _density, mean_velocity.x,
                                            mean_velocity.y);
        }
        density[i] = _density;
        vel_x[i] = mean_velocity.x;
        vel_y[i] = mean_velocity.y;
      }
    }
  }

  // Collision step: for each direction d, relax every node's f toward
  // equilibrium (BGK). Looping direction-outer / node-inner (instead of the
  // old node-outer / direction-inner order) is what makes this
  // vectorizable: each inner iteration reads/writes contiguous floats from
  // old_f[d], density, vel_x, vel_y.
  void apply_collision_step() {
    for (int d = 0; d < 9; d++) {
      float *f = old_f[d].data();
      for (std::size_t i = 0; i < N; i++) {
        const float feq =
            compute_equilibrium(d, density[i], vel_x[i], vel_y[i]);
        f[i] = f[i] - (f[i] - feq) / relaxation_time;
      }
    }
  }

  // Here we are assuming perfect elastic container. That is, velocity is
  // reflected symmetrically if it hits the wall.
  void stream() {
    for (int x = 0; x < SIZE; x++) {
      for (int y = 0; y < SIZE; y++) {
        const std::size_t i = idx(x, y);

        for (int d = 0; d < 9; d++) {
          const int dx = static_cast<int>(D2Q9::ex[d]);
          const int dy = static_cast<int>(D2Q9::ey[d]);
          int prev_x, prev_y, prev_dx = dx, prev_dy = dy;

          if (x - dx < 0 || x - dx >= SIZE) {
            prev_x = x;
            prev_dx = -dx;
          } else {
            prev_x = x - dx;
          }

          if (y - dy < 0 || y - dy >= SIZE) {
            prev_y = y;
            prev_dy = -dy;
          } else {
            prev_y = y - dy;
          }

          const int prev_d = D2Q9::dir_index(prev_dx, prev_dy);
          new_f[d][i] = old_f[prev_d][idx(prev_x, prev_y)];
        }

        float _density = 0, mx = 0, my = 0;
        for (int d = 0; d < 9; d++) {
          _density += new_f[d][i];
          mx += D2Q9::ex[d] * new_f[d][i];
          my += D2Q9::ey[d] * new_f[d][i];
        }
        density[i] = _density;
        vel_x[i] = mx / _density;
        vel_y[i] = my / _density;

        min_density = std::min(min_density, _density);
        max_density = std::max(max_density, _density);
        const float speed =
            std::sqrt(vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i]);
        min_speed = std::min(min_speed, speed);
        max_speed = std::max(max_speed, speed);
      }
    }
    old_f = new_f;
  }

  void apply_step() {
    apply_collision_step();
    stream();
  }
};

int main() {
  const float viscosity = 0.001f;
  const int metres_per_lattice = 1;
  const float u_lattice =
      0.08; // must be less than lattice speed of sound ~0.577
  const float lattice_time_unit = metres_per_lattice / u_lattice;
  // BGK is only stable for tau > 0.5; the viscosity-derived value lands at
  // 0.50024, right on the boundary, so set it directly for now
  const float relaxation_time = 0.6f;
  std::unique_ptr<Lattice> lattice = std::make_unique<Lattice>(relaxation_time);
  Lattice &simulation = *lattice;

  // Using Raylib to simulate Fluid Flow
  const int latticeSize = 2;
  InitWindow(SIZE * latticeSize, SIZE * latticeSize, "Fluuuidsss");
  SetTargetFPS(100);

  while (!WindowShouldClose()) {
    simulation.apply_step();

    const float density_span =
        std::max(simulation.max_density - simulation.min_density, 1e-4f);

    // Drawing Buffer
    BeginDrawing();
    ClearBackground(WHITE);
    for (int x = 0; x < SIZE; x++) {
      for (int y = 0; y < SIZE; y++) {
        const std::size_t i = Lattice::idx(x, y);
        float density = simulation.density[i];
        float scaled_density =
            (density - simulation.min_density) / density_span;

        // If the solver has diverged, clamp rather than feed NaN to
        // ColorFromHSV
        if (!std::isfinite(scaled_density))
          scaled_density = 0.0f;

        Color c = ColorFromHSV(240.0f * (1.0f - scaled_density), 1.0f, 1.0f);
        DrawRectangle(x * latticeSize, y * latticeSize, latticeSize,
                      latticeSize, c);
      }
    }
    DrawFPS(10, 10);
    EndDrawing();
  }

  return 0;
}
