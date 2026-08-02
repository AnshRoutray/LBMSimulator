#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <memory>
#include <raylib.h>
#include <stdint.h>
#include <unistd.h>

constexpr std::size_t SIZE = 300;
constexpr std::size_t DIMENSIONS = 2;
constexpr float cs2 = 1.0 / 3.0;

namespace D2Q9 {
// Gauss Hermite weights_[x+1][y+1] for D2Q9
static constexpr std::array<std::array<double, 3>, 3> weights_ = {
    {{{1.0 / 36, 1.0 / 9, 1.0 / 36}},
     {{1.0 / 9, 4.0 / 9, 1.0 / 9}},
     {{1.0 / 36, 1.0 / 9, 1.0 / 36}}}};

static constexpr double weight(int dx, int dy) {
  return weights_[dx + 1][dy + 1];
}
}; // namespace D2Q9

struct Vector {
  float x, y;
  Vector(float _x, float _y) : x{_x}, y{_y} {}
  Vector() : x{0}, y{0} {}
  Vector operator+(const Vector v2) const { return Vector(x + v2.x, y + v2.y); }
  Vector operator*(const float multiple) const {
    return Vector(x * multiple, y * multiple);
  }
  Vector operator/(const float divisor) const {
    return Vector(x / divisor, y / divisor);
  }
  void operator+=(const Vector v2) { *this = *this + v2; }
  void operator*=(const float multiple) { *this = *this * multiple; }
  void operator/=(const float divisor) { *this = *this / divisor; }
  void zero() { x = 0, y = 0; }
  float dot(const Vector v2) const { return x * v2.x + y * v2.y; }
};

class LatticeNode {
  std::array<float, (DIMENSIONS + 1) * (DIMENSIONS + 1)> _vDistribution;

public:
  float density;
  Vector mean_velocity;
  // This makes it so that the array can be indexed by -1, 0, and 1
  float &vDistribution(int x, int y) {
    return _vDistribution[(x + 1) * 3 + (y + 1)];
  }

  void compute_mean_velocity() {
    mean_velocity.zero();
    for (int x = -1; x <= 1; x++) {
      for (int y = -1; y <= 1; y++) {
        mean_velocity += (Vector(x, y) * vDistribution(x, y));
      }
    }
    mean_velocity /= density;
  }

  void compute_density() {
    density = 0;
    for (int x = -1; x <= 1; x++) {
      for (int y = -1; y <= 1; y++) {
        density += vDistribution(x, y);
      }
    }
  }

  static float compute_equilibrium_velocity(Vector e, float _density,
                                            Vector mean_velocity) {
    float dot_product = e.dot(mean_velocity);
    float new_vDistrNum = (1 + (dot_product / cs2) +
                           (dot_product * dot_product) / (2 * cs2 * cs2) -
                           mean_velocity.dot(mean_velocity) / (2 * cs2)) *
                          _density * D2Q9::weight(e.x, e.y);
    return new_vDistrNum;
  }

  void apply_collision_step(float relaxation_time) {
    // Computing f^eq
    for (int x = -1; x <= 1; x++) {
      for (int y = -1; y <= 1; y++) {
        // Computing Equilibrium Distribution
        float old_vDistrNum = vDistribution(x, y);
        float new_vDistrNum =
            compute_equilibrium_velocity(Vector(x, y), density, mean_velocity);

        // Applying the Bhatnagar-Gross-Kook model to relax the current
        // state towards Equilibrium
        vDistribution(x, y) =
            old_vDistrNum -
            1 / relaxation_time * (old_vDistrNum - new_vDistrNum);
      }
    }
  }
};

class Lattice {
public:
  // Two copies so that streaming phase is not corrupted, when shifting velocity
  // distributions between lattices
  std::array<std::array<LatticeNode, SIZE>, SIZE> new_lattice;
  std::array<std::array<LatticeNode, SIZE>, SIZE> old_lattice;
  float relaxation_time;
  float min_density;
  float max_density;

  // The peak:ambient ratio is a pressure ratio, and BGK is only valid while the
  // flow it drives stays well below the lattice sound speed (cs = 0.577). 1.3:1
  // peaks at about Mach 0.11; 4:1 runs away past Mach 1 and blows up.
  Lattice(float _relaxation_time, float _peak_density = 0.2f,
          float _base_density = 0.15f, float _blob_radius = SIZE * 0.15f,
          Vector _mean_velocity = Vector(0.0f, 0.0f))
      : relaxation_time{_relaxation_time}, min_density{FLT_MAX},
        max_density(-FLT_MAX) {
    init(_peak_density, _base_density, _blob_radius, _mean_velocity);
  }

  // Starts the fluid at rest as a radially symmetric blob: peak_density at the
  // centre of the domain falling off to ambient_density at the edges. The
  // resulting pressure gradient is what drives the flow on the first steps.
  //
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
        // sharp edge launches a ringing pressure wave on the first step
        const float falloff = std::exp(-(dx * dx + dy * dy) / two_r2);
        const float density =
            base_density + (peak_density - base_density) * falloff;

        LatticeNode &node = old_lattice[x][y];
        for (int vx = -1; vx <= 1; vx++) {
          for (int vy = -1; vy <= 1; vy++) {
            node.vDistribution(vx, vy) =
                LatticeNode::compute_equilibrium_velocity(
                    Vector(vx, vy), density, mean_velocity);
          }
        }
        node.density = density;
        node.mean_velocity = mean_velocity;
      }
    }
  }

  // Here we are assuming perfect elastic container. That is, velocity is
  // reflected symmetrically if it hits the wall
  void stream() {
    min_density = FLT_MAX;
    max_density = -FLT_MAX;

    for (int x = 0; x < SIZE; x++) {
      for (int y = 0; y < SIZE; y++) {
        for (int v_x = -1; v_x <= 1; v_x++) {
          for (int v_y = -1; v_y <= 1; v_y++) {
            int prev_x, prev_y, prev_vx, prev_vy;
            if (x - v_x < 0 || x - v_x >= SIZE) {
              prev_x = x;
              prev_vx = -v_x;
            } else {
              prev_x = x - v_x;
              prev_vx = v_x;
            }

            if (y - v_y < 0 || y - v_y >= SIZE) {
              prev_y = y;
              prev_vy = -v_y;
            } else {
              prev_y = y - v_y;
              prev_vy = v_y;
            }

            new_lattice[x][y].vDistribution(v_x, v_y) =
                old_lattice[prev_x][prev_y].vDistribution(prev_vx, prev_vy);
          }
        }
        new_lattice[x][y].compute_density();
        new_lattice[x][y].compute_mean_velocity();
        min_density = std::min(min_density, new_lattice[x][y].density);
        max_density = std::max(max_density, new_lattice[x][y].density);
      }
    }
    old_lattice = new_lattice;
  }

  void apply_step() {
    for (int x = 0; x < SIZE; x++) {
      for (int y = 0; y < SIZE; y++) {
        old_lattice[x][y].apply_collision_step(relaxation_time);
      }
    }
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
  SetTargetFPS(30);

  while (!WindowShouldClose()) {
    simulation.apply_step();

    const float span =
        std::max(simulation.max_density - simulation.min_density, 1e-4f);

    // Drawing Buffer
    BeginDrawing();
    ClearBackground(WHITE);
    for (int x = 0; x < SIZE; x++) {
      for (int y = 0; y < SIZE; y++) {
        float density = simulation.old_lattice[x][y].density;
        float scaled = (density - simulation.min_density) / span;
        // If the solver has diverged, clamp rather than feed NaN to
        // ColorFromHSV;
        if (!std::isfinite(scaled))
          scaled = 0.0f;
        Color c = ColorFromHSV(240.0f * (1.0f - scaled), 1.0f, 1.0f);
        DrawRectangle(x * latticeSize, y * latticeSize, latticeSize,
                      latticeSize, c);
      }
    }
    DrawFPS(10, 10);
    EndDrawing();
  }

  return 0;
}
