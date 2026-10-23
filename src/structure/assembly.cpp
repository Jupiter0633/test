#include "structure/assembly.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <set>
#include <stdexcept>

namespace hcw::structure {

Model::Model(Mesh mesh, Material mat)
    : mesh_(std::move(mesh)), mat_(mat), fixed_(2 * mesh_.nodes.size(), false) {}

void Model::fix(int node, int comp) {
  if (finalized_) throw std::logic_error("Model::fix after finalize");
  const int d = 2 * node + comp;
  if (!fixed_[d]) {
    fixed_[d] = true;
    fixed_list_.push_back(d);
  }
}

void Model::fix_where(const std::function<bool(double, double)>& pred, int comp) {
  for (int n = 0; n < static_cast<int>(mesh_.nodes.size()); ++n) {
    if (pred(mesh_.nodes[n].r, mesh_.nodes[n].z)) fix(n, comp);
  }
}

void Model::fix_axis() {
  fix_where([](double r, double) { return r == 0.0; }, 0);
}

namespace {

// 절점 그래프의 역 Cuthill–McKee 순서. 연결 성분마다 거의 가장자리인 절점에서 시작한다.
std::vector<int> rcm_order(const Mesh& mesh) {
  const int nn = static_cast<int>(mesh.nodes.size());
  std::vector<std::set<int>> adj(nn);
  for (const auto& e : mesh.elems) {
    for (int a : e) {
      for (int b : e) {
        if (a != b) adj[a].insert(b);
      }
    }
  }
  std::vector<int> order;
  order.reserve(nn);
  std::vector<bool> seen(nn, false);
  auto bfs_levels = [&](int start, std::vector<int>& level) {
    level.assign(nn, -1);
    std::deque<int> q{start};
    level[start] = 0;
    int last = start;
    while (!q.empty()) {
      const int v = q.front();
      q.pop_front();
      last = v;
      for (int w : adj[v]) {
        if (level[w] < 0) {
          level[w] = level[v] + 1;
          q.push_back(w);
        }
      }
    }
    return last;
  };
  std::vector<int> level;
  for (int s0 = 0; s0 < nn; ++s0) {
    if (seen[s0]) continue;
    // 가장 먼 절점을 두 번 따라가 시작점을 정한다
    int start = bfs_levels(s0, level);
    start = bfs_levels(start, level);
    std::vector<int> comp;
    std::deque<int> q{start};
    seen[start] = true;
    while (!q.empty()) {
      const int v = q.front();
      q.pop_front();
      comp.push_back(v);
      std::vector<int> nb;
      for (int w : adj[v]) {
        if (!seen[w]) nb.push_back(w);
      }
      std::sort(nb.begin(), nb.end(), [&](int a, int b) {
        return adj[a].size() != adj[b].size() ? adj[a].size() < adj[b].size() : a < b;
      });
      for (int w : nb) {
        seen[w] = true;
        q.push_back(w);
      }
    }
    order.insert(order.end(), comp.begin(), comp.end());
  }
  std::reverse(order.begin(), order.end());
  return order;
}

}  // namespace

void Model::finalize() {
  if (finalized_) return;
  const int ne = static_cast<int>(mesh_.elems.size());
  ke_.assign(256 * static_cast<std::size_t>(ne), 0.0);
  me_.assign(256 * static_cast<std::size_t>(ne), 0.0);
  for (int k = 0; k < ne; ++k) element_matrices(mesh_, k, mat_, &ke_[256 * k], &me_[256 * k]);

  eq_.assign(num_dofs(), -1);
  dof_of_eq_.clear();
  for (int n : rcm_order(mesh_)) {
    for (int c = 0; c < 2; ++c) {
      const int d = 2 * n + c;
      if (!fixed_[d]) {
        eq_[d] = static_cast<int>(dof_of_eq_.size());
        dof_of_eq_.push_back(d);
      }
    }
  }
  num_free_ = static_cast<int>(dof_of_eq_.size());
  top_.resize(num_free_);
  for (int i = 0; i < num_free_; ++i) top_[i] = i;
  for (const auto& e : mesh_.elems) {
    int lo = num_free_;
    for (int a : e) {
      for (int c = 0; c < 2; ++c) {
        if (eq_[2 * a + c] >= 0) lo = std::min(lo, eq_[2 * a + c]);
      }
    }
    for (int a : e) {
      for (int c = 0; c < 2; ++c) {
        const int q = eq_[2 * a + c];
        if (q >= 0) top_[q] = std::min(top_[q], lo);
      }
    }
  }
  finalized_ = true;
}

void Model::apply(double cK, double cM, const Vec& x, Vec& y) const {
  if (!finalized_) throw std::logic_error("Model::apply before finalize");
  y.assign(num_dofs(), 0.0);
  double xe[16];
  for (int k = 0; k < static_cast<int>(mesh_.elems.size()); ++k) {
    const auto& e = mesh_.elems[k];
    for (int a = 0; a < 8; ++a) {
      xe[2 * a] = x[2 * e[a]];
      xe[2 * a + 1] = x[2 * e[a] + 1];
    }
    const double* K = &ke_[256 * static_cast<std::size_t>(k)];
    const double* M = &me_[256 * static_cast<std::size_t>(k)];
    for (int r = 0; r < 16; ++r) {
      double s = 0.0;
      for (int c = 0; c < 16; ++c) s += (cK * K[16 * r + c] + cM * M[16 * r + c]) * xe[c];
      y[2 * e[r / 2] + r % 2] += s;
    }
  }
}

void Model::apply_K(const Vec& x, Vec& y) const { apply(1.0, 0.0, x, y); }
void Model::apply_M(const Vec& x, Vec& y) const { apply(0.0, 1.0, x, y); }

SkylineMatrix Model::assemble(double cK, double cM) const {
  if (!finalized_) throw std::logic_error("Model::assemble before finalize");
  SkylineMatrix A(top_);
  for (int k = 0; k < static_cast<int>(mesh_.elems.size()); ++k) {
    const auto& e = mesh_.elems[k];
    int q[16];
    for (int a = 0; a < 8; ++a) {
      q[2 * a] = eq_[2 * e[a]];
      q[2 * a + 1] = eq_[2 * e[a] + 1];
    }
    const double* K = &ke_[256 * static_cast<std::size_t>(k)];
    const double* M = &me_[256 * static_cast<std::size_t>(k)];
    for (int r = 0; r < 16; ++r) {
      if (q[r] < 0) continue;
      for (int c = 0; c < 16; ++c) {
        if (q[c] < 0 || q[c] < q[r]) continue;  // 위 삼각만 (i ≤ j)
        A.add(q[r], q[c], cK * K[16 * r + c] + cM * M[16 * r + c]);
      }
    }
  }
  return A;
}

Vec Model::restrict_free(const Vec& full) const {
  Vec f(num_free_);
  for (int i = 0; i < num_free_; ++i) f[i] = full[dof_of_eq_[i]];
  return f;
}

void Model::scatter_free(const Vec& free, Vec& full) const {
  full.resize(num_dofs(), 0.0);
  for (int i = 0; i < num_free_; ++i) full[dof_of_eq_[i]] = free[i];
}

double Model::dot(const Vec& a, const Vec& b) const {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

namespace {
double norm2(const Vec& v) {
  double s = 0.0;
  for (double x : v) s += x * x;
  return std::sqrt(s);
}
}  // namespace

Vec solve_checked(const Model& model, const SkylineMatrix& factored, double cK, double cM, const Vec& rhs,
                  double* rel_residual) {
  Vec x = rhs;
  factored.solve(x);
  // 잔차는 저장된 요소 행렬로 따로 계산한다 (분해와 독립).
  // 얇은 판처럼 조건수가 큰 행렬은 잔차 보정(x += A⁻¹ r)을 최대 3 번 한다.
  Vec full(model.num_dofs(), 0.0), y;
  const double nb = norm2(rhs);
  double res = 0.0;
  for (int pass = 0;; ++pass) {
    model.scatter_free(x, full);
    model.apply(cK, cM, full, y);
    const Vec yf = model.restrict_free(y);
    Vec r(rhs.size());
    for (std::size_t i = 0; i < r.size(); ++i) r[i] = rhs[i] - yf[i];
    res = nb > 0.0 ? norm2(r) / nb : norm2(r);
    if (res <= 1e-14 || pass == 3) break;
    factored.solve(r);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] += r[i];
  }
  if (rel_residual) *rel_residual = res;
  return x;
}

Vec solve_static(const Model& model, const Vec& F, const Vec& u_prescribed, double* rel_residual) {
  Vec up(model.num_dofs(), 0.0);
  for (int d : model.fixed_dofs()) up[d] = u_prescribed[d];
  Vec Kup;
  model.apply_K(up, Kup);
  Vec rhs = model.restrict_free(F);
  const Vec kf = model.restrict_free(Kup);
  for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] -= kf[i];
  SkylineMatrix K = model.assemble(1.0, 0.0);
  K.factor();
  const Vec x = solve_checked(model, K, 1.0, 0.0, rhs, rel_residual);
  model.scatter_free(x, up);
  return up;
}

Vec reactions(const Model& model, const Vec& u, const Vec& F) {
  Vec R;
  model.apply_K(u, R);
  for (std::size_t i = 0; i < R.size(); ++i) R[i] -= F[i];
  return R;
}

}  // namespace hcw::structure
