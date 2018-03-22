/* Copyright 2016 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "soleil_mapper.h"

#include "default_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

static LegionRuntime::Logger::Category log_soleil("soleil");

class SoleilMapper : public DefaultMapper
{
public:
  SoleilMapper(MapperRuntime *rt, Machine machine, Processor local,
                const char *mapper_name,
                std::vector<Processor>* loc_procs_list,
                std::vector<Processor>* toc_procs_list,
                std::vector<Processor>* omp_procs_list,
                std::vector<Processor>* io_procs_list,
                std::vector<Memory>* sysmems_list,
                std::vector<Memory>* fbmems_list,
                std::map<Memory, std::vector<Processor> >* sysmem_local_procs,
                std::map<Memory, std::vector<Processor> >* sysmem_local_io_procs,
                std::map<Memory, std::vector<Processor> >* fbmem_local_procs,
                std::map<Processor, Memory>* proc_sysmems,
                std::map<Processor, Memory>* proc_regmems,
                std::map<Processor, Memory>* proc_fbmems,
                std::map<Processor, Memory>* proc_zcmems);
  virtual void select_task_options(const MapperContext    ctx,
                                   const Task&            task,
                                         TaskOptions&     output);
  virtual void default_policy_rank_processor_kinds(
                                    MapperContext ctx, const Task &task,
                                    std::vector<Processor::Kind> &ranking);
  virtual void default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs);
  virtual Memory default_policy_select_target_memory(MapperContext ctx,
                                    Processor target_proc,
                                    const RegionRequirement &req);
  virtual LogicalRegion default_policy_select_instance_region(
                                     MapperContext ctx, Memory target_memory,
                                     const RegionRequirement &req,
                                     const LayoutConstraintSet &constraints,
                                     bool force_new_instances,
                                     bool meets_constraints);
  virtual void map_task(const MapperContext      ctx,
                        const Task&              task,
                        const MapTaskInput&      input,
                              MapTaskOutput&     output);
  virtual void map_copy(const MapperContext ctx,
                        const Copy &copy,
                        const MapCopyInput &input,
                        MapCopyOutput &output);
protected:
  template<bool IS_SRC>
  void soleil_create_copy_instance(MapperContext ctx, const Copy &copy,
                                   const RegionRequirement &req, unsigned index,
                                   std::vector<PhysicalInstance> &instances);
private:
  bool use_gpu;
  bool use_omp;
  std::vector<Processor>& loc_procs_list;
  std::vector<Processor>& toc_procs_list;
  std::vector<Processor>& omp_procs_list;
  std::vector<Processor>& io_procs_list;
  std::vector<Memory>& sysmems_list;
  std::vector<Memory>& fbmems_list;
  std::map<Memory, std::vector<Processor> >& sysmem_local_procs;
  std::map<Memory, std::vector<Processor> >& sysmem_local_io_procs;
  std::map<Memory, std::vector<Processor> >& fbmem_local_procs;
  std::map<Processor, Memory>& proc_sysmems;
  std::map<Processor, Memory>& proc_regmems;
  std::map<Processor, Memory>& proc_fbmems;
  std::map<Processor, Memory>& proc_zcmems;
  std::vector<TaskSlice> loc_slice_cache;
  std::vector<TaskSlice> toc_slice_cache;
  std::vector<TaskSlice> omp_slice_cache;
  std::vector<Processor::Kind> ranking;
  std::vector<Processor> initial_procs;
};

//--------------------------------------------------------------------------
SoleilMapper::SoleilMapper(MapperRuntime *rt, Machine machine, Processor local,
                             const char *mapper_name,
                             std::vector<Processor>* _loc_procs_list,
                             std::vector<Processor>* _toc_procs_list,
                             std::vector<Processor>* _omp_procs_list,
                             std::vector<Processor>* _io_procs_list,
                             std::vector<Memory>* _sysmems_list,
                             std::vector<Memory>* _fbmems_list,
                             std::map<Memory, std::vector<Processor> >* _sysmem_local_procs,
                             std::map<Memory, std::vector<Processor> >* _sysmem_local_io_procs,
                             std::map<Memory, std::vector<Processor> >* _fbmem_local_procs,
                             std::map<Processor, Memory>* _proc_sysmems,
                             std::map<Processor, Memory>* _proc_regmems,
                             std::map<Processor, Memory>* _proc_fbmems,
                             std::map<Processor, Memory>* _proc_zcmems)
//--------------------------------------------------------------------------
  : DefaultMapper(rt, machine, local, mapper_name),
    use_gpu(false),
    use_omp(false),
    loc_procs_list(*_loc_procs_list),
    toc_procs_list(*_toc_procs_list),
    omp_procs_list(*_omp_procs_list),
    io_procs_list(*_io_procs_list),
    sysmems_list(*_sysmems_list),
    fbmems_list(*_fbmems_list),
    sysmem_local_procs(*_sysmem_local_procs),
    sysmem_local_io_procs(*_sysmem_local_io_procs),
    fbmem_local_procs(*_fbmem_local_procs),
    proc_sysmems(*_proc_sysmems),
    proc_regmems(*_proc_regmems),
    proc_fbmems(*_proc_fbmems),
    proc_zcmems(*_proc_zcmems)
{
  use_gpu = toc_procs_list.size() > 0;
  use_omp = omp_procs_list.size() > 0;
  if (use_gpu)
  {
    ranking.push_back(Processor::TOC_PROC);
    initial_procs.push_back(*toc_procs_list.begin());
  }
  else if (use_omp)
  {
    ranking.push_back(Processor::OMP_PROC);
    initial_procs.push_back(*omp_procs_list.begin());
  }
  ranking.push_back(Processor::LOC_PROC);
  initial_procs.push_back(*loc_procs_list.begin());
}

void SoleilMapper::select_task_options(const MapperContext    ctx,
                                       const Task&            task,
                                             TaskOptions&     output)
{
  output.initial_proc = default_policy_select_initial_processor(ctx, task);
  output.inline_task = false;
  output.stealable = stealing_enabled;
  output.map_locally = true;
}

void SoleilMapper::default_policy_rank_processor_kinds(MapperContext ctx,
                        const Task &task, std::vector<Processor::Kind> &ranking)
{
  const char* task_name = task.get_task_name();
  const char* prefix = "shard_";
  if (strncmp(task_name, prefix, strlen(prefix)) == 0) {
    // Put shard tasks on IO processors.
    ranking.resize(4);
    ranking[0] = Processor::TOC_PROC;
    ranking[1] = Processor::PROC_SET;
    ranking[2] = Processor::IO_PROC;
    ranking[3] = Processor::LOC_PROC;
  } else {
    DefaultMapper::default_policy_rank_processor_kinds(ctx, task, ranking);
  }
}

void SoleilMapper::default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs)
{
  target_procs.push_back(task.target_proc);
}

static bool is_ghost(MapperRuntime *runtime,
                     const MapperContext ctx,
                     LogicalRegion leaf)
{
  // If the region has no parent then it was from a duplicated
  // partition and therefore must be a ghost.
  if (!runtime->has_parent_logical_partition(ctx, leaf)) {
    return true;
  }

  return false;
}

Memory SoleilMapper::default_policy_select_target_memory(MapperContext ctx,
                                                   Processor target_proc,
                                                   const RegionRequirement &req)
{
  Memory target_memory = proc_sysmems[target_proc];
  if (is_ghost(runtime, ctx, req.region)) {
    std::map<Processor, Memory>::iterator finder = proc_regmems.find(target_proc);
    if (finder != proc_regmems.end()) target_memory = finder->second;
  }
  return target_memory;
}

LogicalRegion SoleilMapper::default_policy_select_instance_region(
                             MapperContext ctx, Memory target_memory,
                             const RegionRequirement &req,
                             const LayoutConstraintSet &layout_constraints,
                             bool force_new_instances,
                             bool meets_constraints)
{
  return req.region;
}

void SoleilMapper::map_task(const MapperContext      ctx,
                            const Task&              task,
                            const MapTaskInput&      input,
                                  MapTaskOutput&     output)
{
  if (task.parent_task != NULL && task.parent_task->must_epoch_task) {
    Processor::Kind target_kind = task.target_proc.kind();
    // Get the variant that we are going to use to map this task
    VariantInfo chosen = default_find_preferred_variant(task, ctx,
                                                        true/*needs tight bound*/, true/*cache*/, target_kind);
    output.chosen_variant = chosen.variant;
    // TODO: some criticality analysis to assign priorities
    output.task_priority = 0;
    output.postmap_task = false;
    // Figure out our target processors
    output.target_procs.push_back(task.target_proc);

    for (unsigned idx = 0; idx < task.regions.size(); idx++) {
      const RegionRequirement &req = task.regions[idx];

      // Skip any empty regions
      if ((req.privilege == NO_ACCESS) || (req.privilege_fields.empty()))
        continue;

      assert(input.valid_instances[idx].size() > 0);
      output.chosen_instances[idx] = input.valid_instances[idx];
      bool ok = runtime->acquire_and_filter_instances(ctx, output.chosen_instances);
      if (!ok) {
        log_soleil.error("failed to acquire instances");
        assert(false);
      }
    }
    return;
  }

  DefaultMapper::map_task(ctx, task, input, output);
}

void SoleilMapper::map_copy(const MapperContext ctx,
                             const Copy &copy,
                             const MapCopyInput &input,
                             MapCopyOutput &output)
{
  log_soleil.spew("Soleil mapper map_copy");
  for (unsigned idx = 0; idx < copy.src_requirements.size(); idx++)
  {
    // Always use a virtual instance for the source.
    output.src_instances[idx].clear();
    output.src_instances[idx].push_back(
      PhysicalInstance::get_virtual_instance());

    // Place the destination instance on the remote node.
    output.dst_instances[idx].clear();
    if (!copy.dst_requirements[idx].is_restricted()) {
      // Call a customized method to create an instance on the desired node.
      soleil_create_copy_instance<false/*is src*/>(ctx, copy,
        copy.dst_requirements[idx], idx, output.dst_instances[idx]);
    } else {
      // If it's restricted, just take the instance. This will only
      // happen inside the shard task.
      output.dst_instances[idx] = input.dst_instances[idx];
      if (!output.dst_instances[idx].empty())
        runtime->acquire_and_filter_instances(ctx,
                                output.dst_instances[idx]);
    }
  }
}

template<bool IS_SRC>
void SoleilMapper::soleil_create_copy_instance(MapperContext ctx,
                     const Copy &copy, const RegionRequirement &req,
                     unsigned idx, std::vector<PhysicalInstance> &instances)
{
  // This method is identical to the default version except that it
  // chooses an intelligent memory based on the destination of the
  // copy.

  // See if we have all the fields covered
  std::set<FieldID> missing_fields = req.privilege_fields;
  for (std::vector<PhysicalInstance>::const_iterator it =
        instances.begin(); it != instances.end(); it++)
  {
    it->remove_space_fields(missing_fields);
    if (missing_fields.empty())
      break;
  }
  if (missing_fields.empty())
    return;
  // If we still have fields, we need to make an instance
  // We clearly need to take a guess, let's see if we can find
  // one of our instances to use.


  const LogicalRegion& region = copy.src_requirements[idx].region;
  IndexPartition ip = runtime->get_parent_index_partition(ctx, region.get_index_space());
  Domain domain = runtime->get_index_partition_color_space(ctx, ip);
  DomainPoint point =
    runtime->get_logical_region_color_point(ctx, region);
  coord_t size_x = domain.rect_data[3] - domain.rect_data[0] + 1;
  coord_t size_y = domain.rect_data[4] - domain.rect_data[1] + 1;
  Color color = point.point_data[0] +
                point.point_data[1] * size_x +
                point.point_data[2] * size_x * size_y;
  Processor proc = Processor::NO_PROC;
  if (use_gpu) {
    proc = toc_procs_list[color % toc_procs_list.size()];
  } else if (use_omp) {
    proc = omp_procs_list[color % omp_procs_list.size()];
  } else {
    proc = loc_procs_list[color % loc_procs_list.size()];
  }
  Memory target_memory = default_policy_select_target_memory(ctx, proc, req);

  bool force_new_instances = false;
  LayoutConstraintSet creation_constraints;
  default_policy_select_constraints(ctx, creation_constraints, target_memory, req);
  creation_constraints.add_constraint(
      FieldConstraint(missing_fields,
                      false/*contig*/, false/*inorder*/));
  instances.resize(instances.size() + 1);
  if (!default_make_instance(ctx, target_memory,
        creation_constraints, instances.back(),
        COPY_MAPPING, force_new_instances, true/*meets*/, req))
  {
    // If we failed to make it that is bad
    log_soleil.error("Soleil mapper failed allocation for "
                     "%s region requirement %d of explicit "
                     "region-to-region copy operation in task %s "
                     "(ID %lld) in memory " IDFMT " for processor "
                     IDFMT ". This means the working set of your "
                     "application is too big for the allotted "
                     "capacity of the given memory under the default "
                     "mapper's mapping scheme. You have three "
                     "choices: ask Realm to allocate more memory, "
                     "write a custom mapper to better manage working "
                     "sets, or find a bigger machine. Good luck!",
                     IS_SRC ? "source" : "destination", idx,
                     copy.parent_task->get_task_name(),
                     copy.parent_task->get_unique_id(),
		                 target_memory.id,
		                 copy.parent_task->current_proc.id);
    assert(false);
  }
}

static void create_mappers(Machine machine,
                           HighLevelRuntime *runtime,
                           const std::set<Processor> &local_procs)
{
  std::vector<Processor>* loc_procs_list = new std::vector<Processor>();
  std::vector<Processor>* toc_procs_list = new std::vector<Processor>();
  std::vector<Processor>* omp_procs_list = new std::vector<Processor>();
  std::vector<Processor>* io_procs_list = new std::vector<Processor>();
  std::vector<Memory>* sysmems_list = new std::vector<Memory>();
  std::vector<Memory>* fbmems_list = new std::vector<Memory>();
  std::map<Memory, std::vector<Processor> >* sysmem_local_procs =
    new std::map<Memory, std::vector<Processor> >();
  std::map<Memory, std::vector<Processor> >* sysmem_local_io_procs =
    new std::map<Memory, std::vector<Processor> >();
  std::map<Memory, std::vector<Processor> >* fbmem_local_procs =
    new std::map<Memory, std::vector<Processor> >();
  std::map<Processor, Memory>* proc_sysmems = new std::map<Processor, Memory>();
  std::map<Processor, Memory>* proc_regmems = new std::map<Processor, Memory>();
  std::map<Processor, Memory>* proc_fbmems = new std::map<Processor, Memory>();
  std::map<Processor, unsigned> proc_fbmems_affinity;
  std::map<Processor, Memory>* proc_zcmems = new std::map<Processor, Memory>();

  std::vector<Machine::ProcessorMemoryAffinity> proc_mem_affinities;
  machine.get_proc_mem_affinity(proc_mem_affinities);

  for (unsigned idx = 0; idx < proc_mem_affinities.size(); ++idx) {
    Machine::ProcessorMemoryAffinity& affinity = proc_mem_affinities[idx];
    if (affinity.p.kind() == Processor::LOC_PROC ||
        affinity.p.kind() == Processor::IO_PROC ||
        affinity.p.kind() == Processor::OMP_PROC) {
      if (affinity.m.kind() == Memory::SYSTEM_MEM) {
        (*proc_sysmems)[affinity.p] = affinity.m;
      }
      else if (affinity.m.kind() == Memory::REGDMA_MEM) {
        (*proc_regmems)[affinity.p] = affinity.m;
      }
    }
    else if (affinity.p.kind() == Processor::TOC_PROC) {
      if (affinity.m.kind() == Memory::GPU_FB_MEM) {
        std::map<Processor, unsigned>::iterator finder =
          proc_fbmems_affinity.find(affinity.p);
        if (finder == proc_fbmems_affinity.end() ||
            finder->second > affinity.latency)
        {
          (*proc_fbmems)[affinity.p] = affinity.m;
          proc_fbmems_affinity[affinity.p] = affinity.latency;
        }
      }
      else if (affinity.m.kind() == Memory::Z_COPY_MEM) {
        (*proc_zcmems)[affinity.p] = affinity.m;
      }
    }
  }

  for (std::map<Processor, Memory>::iterator it = proc_sysmems->begin();
       it != proc_sysmems->end(); ++it) {
    if (it->first.kind() == Processor::LOC_PROC) {
      loc_procs_list->push_back(it->first);
      (*sysmem_local_procs)[it->second].push_back(it->first);
    }
    else if (it->first.kind() == Processor::IO_PROC) {
      (*sysmem_local_io_procs)[it->second].push_back(it->first);
      io_procs_list->push_back(it->first);
    }
    else if (it->first.kind() == Processor::OMP_PROC) {
      omp_procs_list->push_back(it->first);
    }
  }

  for (std::map<Memory, std::vector<Processor> >::iterator it =
        sysmem_local_procs->begin(); it != sysmem_local_procs->end(); ++it)
    sysmems_list->push_back(it->first);

  for (std::map<Processor, Memory>::iterator it = proc_fbmems->begin();
       it != proc_fbmems->end(); ++it) {
    if (it->first.kind() == Processor::TOC_PROC) {
      toc_procs_list->push_back(it->first);
      (*fbmem_local_procs)[it->second].push_back(it->first);
    }
  }

  for (std::map<Memory, std::vector<Processor> >::iterator it =
        fbmem_local_procs->begin(); it != fbmem_local_procs->end(); ++it)
    fbmems_list->push_back(it->first);

  for (std::set<Processor>::const_iterator it = local_procs.begin();
        it != local_procs.end(); it++)
  {
    SoleilMapper* mapper = new SoleilMapper(runtime->get_mapper_runtime(),
                                            machine, *it, "soleil_mapper",
                                            loc_procs_list,
                                            toc_procs_list,
                                            omp_procs_list,
                                            io_procs_list,
                                            sysmems_list,
                                            fbmems_list,
                                            sysmem_local_procs,
                                            sysmem_local_io_procs,
                                            fbmem_local_procs,
                                            proc_sysmems,
                                            proc_regmems,
                                            proc_fbmems,
                                            proc_zcmems);
    runtime->replace_default_mapper(mapper, *it);
  }
}

void register_mappers()
{
  HighLevelRuntime::set_registration_callback(create_mappers);
}
