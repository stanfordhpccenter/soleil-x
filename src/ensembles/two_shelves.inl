#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <assert.h>

//struct Config
//{
//  Config(int _x, int _y, int _z, int _t)
//    : x(_x), y(_y), z(_z), t(_t)
//  {}
//  int x;
//  int y;
//  int z;
//  int t;
//};

template<typename T,
         double (*COST_FN)(const T&, unsigned)>
unsigned min_cores_for_target(const T &config, unsigned total_num_cores,
                              double target)
{
  for (unsigned cores = 1; cores <= total_num_cores; ++cores)
    if (COST_FN(config, cores) <= target) return cores;
  return total_num_cores * 2;
}

enum Shelf
{
  S1 = 0,
  S2 = 1
};

struct TaskMapping
{
  TaskMapping()
    : main_idx(0), start_idx(0), end_idx(0)
  {}
  TaskMapping(unsigned m, unsigned s, unsigned e)
    : main_idx(m), start_idx(s), end_idx(e)
  {}
  TaskMapping(const TaskMapping &rhs)
    : main_idx(rhs.main_idx), start_idx(rhs.start_idx), end_idx(rhs.end_idx)
  {}
  ~TaskMapping()
  {}
  unsigned main_idx;
  unsigned start_idx;
  unsigned end_idx;
};

inline
TaskMapping compute_mapping(std::map<unsigned, unsigned> &main_assigned,
                        std::vector<unsigned> &num_analysis_tasks,
                        unsigned cores_per_node, unsigned num_nodes,
                        unsigned start_idx, unsigned end_idx,
                        unsigned num_cores = 0)
{
  if (num_cores == 0) num_cores = end_idx - start_idx + 1;
  unsigned start_node_id = start_idx / cores_per_node;
  unsigned end_node_id = end_idx / cores_per_node;

  unsigned min_num_analysis_tasks = num_analysis_tasks[start_node_id];
  unsigned min_node_id = start_node_id;
  for (unsigned id = start_node_id; id <= end_node_id; ++id)
    if (min_num_analysis_tasks > num_analysis_tasks[id])
    {
      min_num_analysis_tasks = num_analysis_tasks[id];
      min_node_id = id;
    }

  unsigned search_main_start =
    std::max(start_idx, min_node_id * cores_per_node);
  unsigned search_main_end =
    std::min(end_idx, (min_node_id + 1) * cores_per_node - 1);
  num_analysis_tasks[min_node_id] += num_cores;

  unsigned min_num_mains = main_assigned.size();
  unsigned min_num_mains_idx = search_main_start;
  for (unsigned idx = search_main_start; idx <= search_main_end; ++idx)
  {
    std::map<unsigned, unsigned>::iterator finder = main_assigned.find(idx);
    if (finder == main_assigned.end())
    {
      min_num_mains_idx = idx;
      main_assigned[idx] = 0;
      break;
    }
    else if (finder->second < min_num_mains)
    {
      min_num_mains_idx = idx;
      min_num_mains = finder->second;
    }
  }
  ++main_assigned[min_num_mains_idx];

  assert(start_idx <= min_num_mains_idx && min_num_mains_idx <= end_idx);

  return TaskMapping(min_num_mains_idx, start_idx, end_idx);
}

template<typename T,
         double (*COST_FN)(const T&, unsigned)>
void two_shelves(const std::vector<T> &configs,
                 std::vector<TaskMapping> &mappings,
                 unsigned cores_per_node, unsigned num_nodes)
{
  unsigned total_num_cores = cores_per_node * num_nodes;
  double min_cost = 0;
  double max_cost = 0;
  for (typename std::vector<T>::const_iterator it = configs.begin();
       it != configs.end(); ++it)
  {
    min_cost += COST_FN(*it, 1);
    max_cost += COST_FN(*it, total_num_cores);
  }
  min_cost /= total_num_cores;

  double target = (min_cost + max_cost) / 2;

  for (unsigned idx = 0; idx < configs.size(); ++idx)
    mappings.push_back(TaskMapping(idx, 0, total_num_cores - 1));

  // Implementation of two-shelf scheduling algorithm as per Mounie et al.
  while (std::fabs(min_cost - max_cost) / max_cost > 1e-5)
  {
    std::vector<unsigned> small_tasks;
    std::map<unsigned, unsigned> cores_small_tasks;
    std::vector<unsigned> large_tasks;
    std::map<unsigned, unsigned> main_assigned;
    std::vector<unsigned> num_analysis_tasks;
    num_analysis_tasks.resize(num_nodes);
    for (unsigned idx = 0; idx < num_nodes; ++idx)
      num_analysis_tasks[idx] = 0;
    for (unsigned idx = 0; idx < configs.size(); ++idx)
      //if (COST_FN(configs[idx], 1) * 2 <= target)
      //{
      //  small_tasks.push_back(idx);
      //  cores_small_tasks[idx] = 1;
      //}
      //else
        large_tasks.push_back(idx);

    double W[large_tasks.size() + 1][total_num_cores + 1];
    int S[large_tasks.size() + 1][total_num_cores + 1];
    unsigned C[large_tasks.size() + 1][total_num_cores + 1];
    bool is_max[large_tasks.size() + 1][total_num_cores + 1];
    W[0][0] = 0;
    S[0][0] = S1;
    C[0][0] = 0;
    for (unsigned q = 1; q <= total_num_cores; ++q)
    {
      W[0][q] = 0;
      S[0][q] = S1;
      C[0][q] = 0;
      is_max[0][q] = false;
    }
    for (unsigned j = 1; j <= large_tasks.size(); ++j)
    {
      W[j][0] = std::numeric_limits<double>::max();
      S[j][0] = S1;
      C[j][0] = total_num_cores * 2;
      is_max[j][0] = true;
    }
    is_max[0][0] = false;
    for (unsigned j = 1; j <= large_tasks.size(); ++j)
    {
      const Config &config = configs[large_tasks[j - 1]];
      for (unsigned q = 1; q <= total_num_cores; ++q)
      {
        unsigned c1 =
          min_cores_for_target<T, COST_FN>(config, total_num_cores, target);
        double w1 = std::numeric_limits<double>::max();
        bool m1 = true;
        if (q >= c1 && !is_max[j - 1][q - c1])
        {
          w1 = W[j - 1][q - c1] + c1 * COST_FN(config, c1);
          m1 = false;
        }

        unsigned c2 =
          min_cores_for_target<T, COST_FN>(config, total_num_cores, target / 2);
        double w2 = std::numeric_limits<double>::max();
        bool m2 = true;
        if (c2 <= total_num_cores && !is_max[j - 1][q])
        {
          w2 = W[j - 1][q] + c2 * COST_FN(config, c2);
          m2 = false;
        }

        if (w1 < w2)
        {
          W[j][q] = w1;
          S[j][q] = S1;
          C[j][q] = c1;
          is_max[j][q] = m1;
        }
        else
        {
          W[j][q] = w2;
          S[j][q] = S2;
          C[j][q] = c2;
          is_max[j][q] = m2;
        }
      }
    }

    std::vector<unsigned> tasks_S0;
    std::map<unsigned, unsigned> cores_S0;
    unsigned num_cores_S0 = 0;
    std::vector<unsigned> tasks_S1;
    std::map<unsigned, unsigned> cores_S1;
    std::vector<unsigned> tasks_S2;
    std::map<unsigned, unsigned> cores_S2;
    double total_cost_S1 = 0.0;
    unsigned num_cores_S1 = 0;
    unsigned num_cores_S2 = 0;
    unsigned q = total_num_cores;
    for (unsigned j = large_tasks.size(); j >= 1; --j)
    {
      if (is_max[j][q])
      {
        total_cost_S1 = std::numeric_limits<double>::max();
        break;
      }
      unsigned c = C[j][q];
      unsigned task_id = large_tasks[j - 1];
      if (S[j][q] == S1)
      {
        if (q < c)
        {
          total_cost_S1 = std::numeric_limits<double>::max();
          break;
        }
        total_cost_S1 += COST_FN(configs[task_id], c);
        num_cores_S1 += c;
        tasks_S1.push_back(task_id);
        cores_S1[task_id] = c;
        q -= c;
      }
      else
      {
        num_cores_S2 += c;
        tasks_S2.push_back(task_id);
        cores_S2[task_id] = c;
      }
    }

    bool found_feasible = false;
    // If tasks in S1 can never fit to the target, no feasible schedule exists
    if (total_cost_S1 > target * total_num_cores)
      found_feasible = false;
    // If tasks in S2 already fit to a given set of cores, the schedule can be
    // feasible.
    else if (num_cores_S2 <= total_num_cores)
    {
      // If small tasks also can be inserted to the empty space in the schedule,
      // we found a feasible schedule.
      if (small_tasks.size() <= 2 * (total_num_cores - num_cores_S1) +
          (total_num_cores - num_cores_S2))
        found_feasible = true;
    }
    else
    {
      while (!found_feasible)
      {
        num_cores_S2 = 0;
        std::set<unsigned> S2_to_S0;
        unsigned remaining_cores = total_num_cores - num_cores_S0;
        for (std::vector<unsigned>::iterator it = tasks_S2.begin();
             it != tasks_S2.end(); ++it)
          if (cores_S2[*it] < remaining_cores)
          {
            unsigned c = cores_S2[*it];
            remaining_cores -= c;
            num_cores_S2 += c;
          }
          else
          {
            const Config &config = configs[*it];
            unsigned c =
              min_cores_for_target<T, COST_FN>(config, total_num_cores, 3 * target / 2);
            S2_to_S0.insert(*it);
            cores_S0[*it] = c;
            num_cores_S0 += c;
          }
        std::vector<unsigned> new_tasks_S2;
        for (std::vector<unsigned>::iterator it = tasks_S2.begin();
             it != tasks_S2.end(); ++it)
        {
          if (S2_to_S0.find(*it) != S2_to_S0.end())
            tasks_S0.push_back(*it);
          else
            new_tasks_S2.push_back(*it);
        }
        tasks_S2.swap(new_tasks_S2);
        unsigned num_cores_two_shelves = std::max(num_cores_S1, num_cores_S2);
        if (num_cores_two_shelves + num_cores_S0 <= total_num_cores)
          found_feasible = true;
        else
        {
          if (num_cores_S0 > total_num_cores) break;
          if (num_cores_two_shelves > 0 && tasks_S2.size() == 0) break;
        }
      }
      // Insert small tasks to the schedule
      if (found_feasible)
      {
        unsigned remaining_cores = total_num_cores -
          (std::max(num_cores_S1, num_cores_S2) + num_cores_S0);
        unsigned num_cores_small = small_tasks.size();
              // Insert small tasks to S1
        if (!((num_cores_small + 1) / 2 + num_cores_S1 <= num_cores_S2 ||
              // Insert small tasks to S2
              num_cores_small + num_cores_S2 <= num_cores_S1 ||
              // Insert small tasks to idle cores
              remaining_cores + num_cores_small <= total_num_cores))
          // If all fail, no feasible schedule exists
          found_feasible = false;
      }
    }

    if (!found_feasible)
    {
      min_cost = target;
      target = (target + max_cost) / 2;
    }
    else
    {
      mappings.clear();
      mappings.resize(configs.size());
      unsigned cores = 0;
      for (std::vector<unsigned>::iterator it = tasks_S0.begin();
          it != tasks_S0.end(); ++it)
      {
        unsigned assigned_cores = cores_S0[*it];
        mappings[*it] = compute_mapping(main_assigned, num_analysis_tasks,
                                        cores_per_node, num_nodes,
                                        cores, cores + assigned_cores - 1);
        cores += assigned_cores;
      }
      assert(cores == num_cores_S0);
      for (std::vector<unsigned>::iterator it = tasks_S1.begin();
          it != tasks_S1.end(); ++it)
      {
        unsigned assigned_cores = cores_S1[*it];
        mappings[*it] = compute_mapping(main_assigned, num_analysis_tasks,
                                        cores_per_node, num_nodes,
                                        cores, cores + assigned_cores - 1);
        cores += assigned_cores;
      }
      assert(cores - num_cores_S0 == num_cores_S1);

      //if (small_tasks.size() == 0)
      {
        std::vector<unsigned>::iterator it = tasks_S2.begin();
        for (; it != tasks_S2.end(); ++it)
        {
          unsigned assigned_cores = cores_S2[*it];
          if (cores + assigned_cores - 1 >= total_num_cores) break;
          mappings[*it] = compute_mapping(main_assigned, num_analysis_tasks,
                                          cores_per_node, num_nodes,
                                          cores, cores + assigned_cores - 1);
          cores += assigned_cores;
        }
        cores = num_cores_S0;
        for (; it != tasks_S2.end(); ++it)
        {
          unsigned assigned_cores = cores_S2[*it];
          mappings[*it] = compute_mapping(main_assigned, num_analysis_tasks,
                                          cores_per_node, num_nodes,
                                          cores, cores + assigned_cores - 1);
          cores += assigned_cores;
        }
      }
      //else
      //{
      //  for (std::vector<unsigned>::iterator it = tasks_S2.begin();
      //       it != tasks_S2.end(); ++it)
      //  {
      //    unsigned assigned_cores = cores_S2[*it];
      //    mappings[*it] = compute_mapping(main_assigned, num_analysis_tasks,
      //                                    cores_per_node, num_nodes,
      //                                    cores, cores + assigned_cores - 1);
      //    cores += assigned_cores;
      //  }
      //  assert(cores - num_cores_S0 == num_cores_S2);
      //}
      //unsigned num_cores_two_shelves = std::max(num_cores_S1, num_cores_S2);
      //cores = num_cores_S0 + num_cores_two_shelves;
      //std::vector<unsigned>::iterator it = small_tasks.begin();
      //for (; it != small_tasks.end(); ++it)
      //{
      //  if (cores >= total_num_cores) break;
      //  mappings[*it] = TaskMapping(cores, cores, cores);
      //  ++cores;
      //}
      //unsigned max_cores = num_cores_S0 + num_cores_two_shelves;
      //for (; it != small_tasks.end(); ++it)
      //{
      //  TaskMapping mapping =
      //    compute_mapping(main_assigned, num_analysis_tasks, cores_per_node,
      //                    num_nodes,
      //                    num_cores_S0 + std::min(num_cores_S1, num_cores_S2),
      //                    max_cores - 1, 1);
      //  mappings[*it] = TaskMapping(mapping.main_idx, mapping.main_idx,
      //                          mapping.main_idx);
      //}
      //assert(it == small_tasks.end());

      max_cost = target;
      target = (target + min_cost) / 2;
    }
  }

#if 0
  for (unsigned idx = 0; idx < mappings.size(); ++idx)
  {
    TaskMapping &mapping = mappings[idx];
    //if (mapping.end_idx - mapping.start_idx + 1 <= cores_per_node)
    //  continue;
    unsigned start_node_id = mapping.start_idx / cores_per_node;
    unsigned end_node_id = mapping.end_idx / cores_per_node;
    unsigned left_spill =
      (start_node_id + 1) * cores_per_node - mapping.start_idx;
    unsigned right_spill =
      mapping.end_idx - end_node_id * cores_per_node;
    fprintf(stderr, "start node %u end node %u left spill %u right spill %u\n",
        start_node_id, end_node_id,left_spill, right_spill);
    //if (left_spill * 8 <= cores_per_node)
    //  mapping.start_idx += left_spill;
    //if (right_spill * 8 <= cores_per_node)
    //  mapping.end_idx -= right_spill;
  }
#endif

  //fprintf(stderr, "Final schedule:\n");
  //for (unsigned idx = 0; idx < mappings.size(); ++idx)
  //  fprintf(stderr, "T%u: %u -- %u (main: %u, #cores: %u)\n",
  //      idx, mappings[idx].start_idx, mappings[idx].end_idx,
  //      mappings[idx].main_idx,
  //      mappings[idx].end_idx - mappings[idx].start_idx + 1);
}

//int main()
//{
//  std::vector<Config> configs;
//  std::vector<TaskMapping> mappings;
//  int cores_per_node = 8;
//  int num_nodes = 4;
//  configs.push_back(Config( 64,  64,  64, 20));
//  configs.push_back(Config( 64,  64,  64, 20));
//  configs.push_back(Config(128,  64,  64, 20));
//  configs.push_back(Config(128, 128, 128, 10));
//  configs.push_back(Config(256, 128, 128, 20));
//  configs.push_back(Config(256, 256, 256, 10));
//
//  two_shelf<Config, cost_fn>(configs, mappings, cores_per_node, num_nodes);
//
//  return 0;
//}
