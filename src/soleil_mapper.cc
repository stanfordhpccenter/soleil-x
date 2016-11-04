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
                std::vector<Processor>* procs_list,
                std::vector<Memory>* sysmems_list,
                std::map<Memory, std::vector<Processor> >* sysmem_local_procs,
                std::map<Processor, Memory>* proc_sysmems,
                std::map<Processor, Memory>* proc_fbmems,
                std::map<Processor, Memory>* proc_zcmems);
  virtual Processor default_policy_select_initial_processor(
                                    MapperContext ctx, const Task &task);
  virtual void default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs);
  virtual void map_task(const MapperContext      ctx,
                        const Task&              task,
                        const MapTaskInput&      input,
                              MapTaskOutput&     output);
  virtual void map_copy(const MapperContext ctx,
                        const Copy &copy,
                        const MapCopyInput &input,
                        MapCopyOutput &output);
protected:
  bool soleil_create_custom_instances(MapperContext ctx,
                          Processor target, Memory target_memory,
                          const RegionRequirement &req, unsigned index,
                          std::set<FieldID> &needed_fields, // will destroy
                          const TaskLayoutConstraintSet &layout_constraints,
                          bool needs_field_constraint_check,
                          std::vector<PhysicalInstance> &instances);
  template<bool IS_SRC>
  void soleil_create_copy_instance(MapperContext ctx, const Copy &copy,
                                   const RegionRequirement &req, unsigned index,
                                   std::vector<PhysicalInstance> &instances);
private:
  std::vector<Processor>& procs_list;
  std::vector<Memory>& sysmems_list;
  std::map<Memory, std::vector<Processor> >& sysmem_local_procs;
  std::map<Processor, Memory>& proc_sysmems;
  std::map<Processor, Memory>& proc_fbmems;
  std::map<Processor, Memory>& proc_zcmems;
};

//--------------------------------------------------------------------------
SoleilMapper::SoleilMapper(MapperRuntime *rt, Machine machine, Processor local,
                             const char *mapper_name,
                             std::vector<Processor>* _procs_list,
                             std::vector<Memory>* _sysmems_list,
                             std::map<Memory, std::vector<Processor> >* _sysmem_local_procs,
                             std::map<Processor, Memory>* _proc_sysmems,
                             std::map<Processor, Memory>* _proc_fbmems,
                             std::map<Processor, Memory>* _proc_zcmems)
//--------------------------------------------------------------------------
  : DefaultMapper(rt, machine, local, mapper_name),
    procs_list(*_procs_list),
    sysmems_list(*_sysmems_list),
    sysmem_local_procs(*_sysmem_local_procs),
    proc_sysmems(*_proc_sysmems),
    proc_fbmems(*_proc_fbmems),
    proc_zcmems(*_proc_zcmems)
{
}

//--------------------------------------------------------------------------
Processor SoleilMapper::default_policy_select_initial_processor(
                                    MapperContext ctx, const Task &task)
//--------------------------------------------------------------------------
{
  return DefaultMapper::default_policy_select_initial_processor(ctx, task);
}

//--------------------------------------------------------------------------
void SoleilMapper::default_policy_select_target_processors(
                                    MapperContext ctx,
                                    const Task &task,
                                    std::vector<Processor> &target_procs)
//--------------------------------------------------------------------------
{
  target_procs.push_back(task.target_proc);
}

//--------------------------------------------------------------------------
void SoleilMapper::map_task(const MapperContext      ctx,
                            const Task&              task,
                            const MapTaskInput&      input,
                                  MapTaskOutput&     output)
//--------------------------------------------------------------------------
{
  Processor::Kind target_kind = task.target_proc.kind();
  // Get the variant that we are going to use to map this task
  VariantInfo chosen = default_find_preferred_variant(task, ctx,
                    true/*needs tight bound*/, true/*cache*/, target_kind);
  output.chosen_variant = chosen.variant;
  // TODO: some criticality analysis to assign priorities
  output.task_priority = 0;
  output.postmap_task = false;
  // Figure out our target processors
  default_policy_select_target_processors(ctx, task, output.target_procs);

  // See if we have an inner variant, if we do virtually map all the regions
  // We don't even both caching these since they are so simple
  if (chosen.is_inner)
  {
    // Check to see if we have any relaxed coherence modes in which
    // case we can no longer do virtual mappings so we'll fall through
    bool has_relaxed_coherence = false;
    for (unsigned idx = 0; idx < task.regions.size(); idx++)
    {
      if (task.regions[idx].prop != EXCLUSIVE)
      {
        has_relaxed_coherence = true;
        break;
      }
    }
    if (!has_relaxed_coherence)
    {
      std::vector<unsigned> reduction_indexes;
      for (unsigned idx = 0; idx < task.regions.size(); idx++)
      {
        // As long as this isn't a reduction-only region requirement
        // we will do a virtual mapping, for reduction-only instances
        // we will actually make a physical instance because the runtime
        // doesn't allow virtual mappings for reduction-only privileges
        if (task.regions[idx].privilege == REDUCE)
          reduction_indexes.push_back(idx);
        else
          output.chosen_instances[idx].push_back(
              PhysicalInstance::get_virtual_instance());
      }
      if (!reduction_indexes.empty())
      {
        const TaskLayoutConstraintSet &layout_constraints =
            runtime->find_task_layout_constraints(ctx,
                                  task.task_id, output.chosen_variant);
        Memory target_memory = default_policy_select_target_memory(ctx,
                                                     task.target_proc);
        for (std::vector<unsigned>::const_iterator it =
              reduction_indexes.begin(); it !=
              reduction_indexes.end(); it++)
        {
          std::set<FieldID> copy = task.regions[*it].privilege_fields;
          if (!soleil_create_custom_instances(ctx, task.target_proc,
              target_memory, task.regions[*it], *it, copy,
              layout_constraints, false/*needs constraint check*/,
              output.chosen_instances[*it]))
          {
            default_report_failed_instance_creation(task, *it,
                                        task.target_proc, target_memory);
          }
        }
      }
      return;
    }
  }
  // First, let's see if we've cached a result of this task mapping
  const unsigned long long task_hash = compute_task_hash(task);
  std::pair<TaskID,Processor> cache_key(task.task_id, task.target_proc);
  std::map<std::pair<TaskID,Processor>,
           std::list<CachedTaskMapping> >::const_iterator
    finder = cached_task_mappings.find(cache_key);
  // This flag says whether we need to recheck the field constraints,
  // possibly because a new field was allocated in a region, so our old
  // cached physical instance(s) is(are) no longer valid
  bool needs_field_constraint_check = false;
  //Memory target_memory = default_policy_select_target_memory(ctx,
  //                                                   task.target_proc);
  Memory target_memory;
  if (task.target_proc.kind() == Processor::LOC_PROC)
    target_memory = proc_sysmems[task.target_proc];
  else if (task.target_proc.kind() == Processor::TOC_PROC)
    target_memory = proc_fbmems[task.target_proc];
  else
    assert(false);

  if (finder != cached_task_mappings.end())
  {
    bool found = false;
    bool has_reductions = false;
    // Iterate through and see if we can find one with our variant and hash
    for (std::list<CachedTaskMapping>::const_iterator it =
          finder->second.begin(); it != finder->second.end(); it++)
    {
      if ((it->variant == output.chosen_variant) &&
          (it->task_hash == task_hash))
      {
        // Have to copy it before we do the external call which
        // might invalidate our iterator
        output.chosen_instances = it->mapping;
        has_reductions = it->has_reductions;
        found = true;
        break;
      }
    }
    if (found)
    {
      // If we have reductions, make those instances now since we
      // never cache the reduction instances
      if (has_reductions)
      {
        const TaskLayoutConstraintSet &layout_constraints =
          runtime->find_task_layout_constraints(ctx,
                              task.task_id, output.chosen_variant);
        for (unsigned idx = 0; idx < task.regions.size(); idx++)
        {
          if (task.regions[idx].privilege == REDUCE)
          {
            std::set<FieldID> copy = task.regions[idx].privilege_fields;
            if (!soleil_create_custom_instances(ctx, task.target_proc,
                target_memory, task.regions[idx], idx, copy,
                layout_constraints, needs_field_constraint_check,
                output.chosen_instances[idx]))
            {
              default_report_failed_instance_creation(task, idx,
                                          task.target_proc, target_memory);
            }
          }
        }
      }
      // See if we can acquire these instances still
      if (runtime->acquire_and_filter_instances(ctx,
                                                 output.chosen_instances))
        return;
      // We need to check the constraints here because we had a
      // prior mapping and it failed, which may be the result
      // of a change in the allocated fields of a field space
      needs_field_constraint_check = true;
      // If some of them were deleted, go back and remove this entry
      // Have to renew our iterators since they might have been
      // invalidated during the 'acquire_and_filter_instances' call
      default_remove_cached_task(ctx, output.chosen_variant,
                    task_hash, cache_key, output.chosen_instances);
    }
  }
  // We didn't find a cached version of the mapping so we need to
  // do a full mapping, we already know what variant we want to use
  // so let's use one of the acceleration functions to figure out
  // which instances still need to be mapped.
  std::vector<std::set<FieldID> > missing_fields(task.regions.size());
  runtime->filter_instances(ctx, task, output.chosen_variant,
                             output.chosen_instances, missing_fields);
  // Track which regions have already been mapped
  std::vector<bool> done_regions(task.regions.size(), false);
  if (!input.premapped_regions.empty())
    for (std::vector<unsigned>::const_iterator it =
          input.premapped_regions.begin(); it !=
          input.premapped_regions.end(); it++)
      done_regions[*it] = true;
  const TaskLayoutConstraintSet &layout_constraints =
    runtime->find_task_layout_constraints(ctx,
                          task.task_id, output.chosen_variant);
  // Now we need to go through and make instances for any of our
  // regions which do not have space for certain fields
  bool has_reductions = false;
  for (unsigned idx = 0; idx < task.regions.size(); idx++)
  {
    if (done_regions[idx])
      continue;
    // Skip any empty regions
    if ((task.regions[idx].privilege == NO_ACCESS) ||
        (task.regions[idx].privilege_fields.empty()) ||
        missing_fields[idx].empty())
      continue;

    if (task.target_proc.kind() == Processor::TOC_PROC &&
        task.regions.size() > 3 && idx >= 3)
      target_memory = proc_zcmems[task.target_proc];

    // See if this is a reduction
    if (task.regions[idx].privilege == REDUCE)
    {
      has_reductions = true;
      if (!soleil_create_custom_instances(ctx, task.target_proc,
              target_memory, task.regions[idx], idx, missing_fields[idx],
              layout_constraints, needs_field_constraint_check,
              output.chosen_instances[idx]))
      {
        default_report_failed_instance_creation(task, idx,
                                    task.target_proc, target_memory);
      }
      continue;
    }
    // Did the application request a virtual mapping for this requirement?
    if ((task.regions[idx].tag & DefaultMapper::VIRTUAL_MAP) != 0)
    {
      PhysicalInstance virt_inst = PhysicalInstance::get_virtual_instance();
      output.chosen_instances[idx].push_back(virt_inst);
      continue;
    }
    // Otherwise make normal instances for the given region
    if (!soleil_create_custom_instances(ctx, task.target_proc,
            target_memory, task.regions[idx], idx, missing_fields[idx],
            layout_constraints, needs_field_constraint_check,
            output.chosen_instances[idx]))
    {
      default_report_failed_instance_creation(task, idx,
                                  task.target_proc, target_memory);
    }
  }
  // Now that we are done, let's cache the result so we can use it later
  std::list<CachedTaskMapping> &map_list = cached_task_mappings[cache_key];
  map_list.push_back(CachedTaskMapping());
  CachedTaskMapping &cached_result = map_list.back();
  cached_result.task_hash = task_hash;
  cached_result.variant = output.chosen_variant;
  cached_result.mapping = output.chosen_instances;
  cached_result.has_reductions = has_reductions;
  // We don't ever save reduction instances in our cache
  if (has_reductions) {
    for (unsigned idx = 0; idx < task.regions.size(); idx++) {
      if (task.regions[idx].privilege != REDUCE)
        continue;
      cached_result.mapping[idx].clear();
    }
  }
}

//--------------------------------------------------------------------------
void SoleilMapper::map_copy(const MapperContext ctx,
                             const Copy &copy,
                             const MapCopyInput &input,
                             MapCopyOutput &output)
//--------------------------------------------------------------------------
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

//--------------------------------------------------------------------------
bool SoleilMapper::soleil_create_custom_instances(MapperContext ctx,
                      Processor target_proc, Memory target_memory,
                      const RegionRequirement &req, unsigned index,
                      std::set<FieldID> &needed_fields,
                      const TaskLayoutConstraintSet &layout_constraints,
                      bool needs_field_constraint_check,
                      std::vector<PhysicalInstance> &instances)
//--------------------------------------------------------------------------
{
  // Before we do anything else figure out our
  // constraints for any instances of this task, then we'll
  // see if these constraints conflict with or are satisfied by
  // any of the other constraints
  bool force_new_instances = false;
  LayoutConstraintID our_layout_id =
   default_policy_select_layout_constraints(ctx, target_memory, req,
           TASK_MAPPING, needs_field_constraint_check, force_new_instances);
  const LayoutConstraintSet &our_constraints =
                runtime->find_layout_constraints(ctx, our_layout_id);
  for (std::multimap<unsigned,LayoutConstraintID>::const_iterator lay_it =
        layout_constraints.layouts.lower_bound(index); lay_it !=
        layout_constraints.layouts.upper_bound(index); lay_it++)
  {
    // Get the constraints
    const LayoutConstraintSet &index_constraints =
              runtime->find_layout_constraints(ctx, lay_it->second);
    std::vector<FieldID> overlaping_fields;
    const std::vector<FieldID> &constraint_fields =
      index_constraints.field_constraint.get_field_set();
    for (unsigned idx = 0; idx < constraint_fields.size(); idx++)
    {
      FieldID fid = constraint_fields[idx];
      std::set<FieldID>::iterator finder = needed_fields.find(fid);
      if (finder != needed_fields.end())
      {
        overlaping_fields.push_back(fid);
        // Remove from the needed fields since we're going to handle it
        needed_fields.erase(finder);
      }
    }
    // If we don't have any overlapping fields, then keep going
    if (overlaping_fields.empty())
      continue;
    // Now figure out how to make an instance
    instances.resize(instances.size()+1);
    // Check to see if these constraints conflict with our constraints
    if (runtime->do_constraints_conflict(ctx,
                                          our_layout_id, lay_it->second))
    {
      // They conflict, so we're just going to make an instance
      // using these constraints
      if (!default_make_instance(ctx, target_memory, index_constraints,
                 instances.back(), TASK_MAPPING, force_new_instances,
                 false/*meets*/, req))
        return false;
    }
    else if (runtime->do_constraints_entail(ctx,
                                             lay_it->second, our_layout_id))
    {
      // These constraints do everything we want to do and maybe more
      // so we can just use them directly
      if (!default_make_instance(ctx, target_memory, index_constraints,
                  instances.back(), TASK_MAPPING, force_new_instances,
                  true/*meets*/, req))
        return false;
    }
    else
    {
      // These constraints don't do as much as we want but don't
      // conflict so make an instance with them and our constraints
      LayoutConstraintSet creation_constraints = index_constraints;
      default_policy_select_constraints(ctx, creation_constraints,
                                        target_memory, req);
      creation_constraints.add_constraint(
          FieldConstraint(overlaping_fields,
            false/*contig*/, false/*inorder*/));
      if (!default_make_instance(ctx, target_memory, creation_constraints,
                     instances.back(), TASK_MAPPING, force_new_instances,
                     true/*meets*/, req))
        return false;
    }
  }
  // If we don't have anymore needed fields, we are done
  if (needed_fields.empty())
    return true;
  // There are no constraints for these fields so we get to do what we want
  instances.resize(instances.size()+1);
  LayoutConstraintSet creation_constraints = our_constraints;
  if (!default_make_instance(ctx, target_memory, creation_constraints,
            instances.back(), TASK_MAPPING, force_new_instances,
            true/*meets*/,  req))
    return false;
  return true;
}

//--------------------------------------------------------------------------
template<bool IS_SRC>
void SoleilMapper::soleil_create_copy_instance(MapperContext ctx,
                     const Copy &copy, const RegionRequirement &req,
                     unsigned idx, std::vector<PhysicalInstance> &instances)
//--------------------------------------------------------------------------
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
  Memory target_memory = default_policy_select_target_memory(ctx,
                           procs_list[color % procs_list.size()]);
  bool force_new_instances = false;
  LayoutConstraintID our_layout_id =
   default_policy_select_layout_constraints(ctx, target_memory,
                                            req, COPY_MAPPING,
                                            true/*needs check*/,
                                            force_new_instances);
  LayoutConstraintSet creation_constraints =
              runtime->find_layout_constraints(ctx, our_layout_id);
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
  std::vector<Processor>* procs_list = new std::vector<Processor>();
  std::vector<Memory>* sysmems_list = new std::vector<Memory>();
  std::map<Memory, std::vector<Processor> >* sysmem_local_procs =
    new std::map<Memory, std::vector<Processor> >();
  std::map<Processor, Memory>* proc_sysmems = new std::map<Processor, Memory>();
  std::map<Processor, Memory>* proc_fbmems = new std::map<Processor, Memory>();
  std::map<Processor, Memory>* proc_zcmems = new std::map<Processor, Memory>();

  std::vector<Machine::ProcessorMemoryAffinity> proc_mem_affinities;
  machine.get_proc_mem_affinity(proc_mem_affinities);

  for (unsigned idx = 0; idx < proc_mem_affinities.size(); ++idx) {
    Machine::ProcessorMemoryAffinity& affinity = proc_mem_affinities[idx];
    if (affinity.p.kind() == Processor::LOC_PROC) {
      if (affinity.m.kind() == Memory::SYSTEM_MEM) {
        (*proc_sysmems)[affinity.p] = affinity.m;
      }
    }
    else if (affinity.p.kind() == Processor::TOC_PROC) {
      if (affinity.m.kind() == Memory::GPU_FB_MEM) {
        (*proc_fbmems)[affinity.p] = affinity.m;
      }
      else if (affinity.m.kind() == Memory::Z_COPY_MEM) {
        (*proc_zcmems)[affinity.p] = affinity.m;
      }
    }
  }

  for (std::map<Processor, Memory>::iterator it = proc_sysmems->begin();
       it != proc_sysmems->end(); ++it) {
    procs_list->push_back(it->first);
    (*sysmem_local_procs)[it->second].push_back(it->first);
  }

  for (std::map<Memory, std::vector<Processor> >::iterator it =
        sysmem_local_procs->begin(); it != sysmem_local_procs->end(); ++it)
    sysmems_list->push_back(it->first);

  for (std::set<Processor>::const_iterator it = local_procs.begin();
        it != local_procs.end(); it++)
  {
    SoleilMapper* mapper = new SoleilMapper(runtime->get_mapper_runtime(),
                                            machine, *it, "soleil_mapper",
                                            procs_list,
                                            sysmems_list,
                                            sysmem_local_procs,
                                            proc_sysmems,
                                            proc_fbmems,
                                            proc_zcmems);
    runtime->replace_default_mapper(mapper, *it);
  }
}

void register_mappers()
{
  HighLevelRuntime::set_registration_callback(create_mappers);
}
