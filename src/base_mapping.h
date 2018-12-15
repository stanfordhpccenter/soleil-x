#ifndef __BASE_MAPPING_H__
#define __BASE_MAPPING_H__

#include <iostream>
#include <stdlib.h>

#include "legion.h"

#include "config_schema.h"

class BaseMapping {
public:
  BaseMapping(const Config& config, unsigned sample_id,
              AddressSpace next_rank, unsigned next_proc_id,
              unsigned prev_proc_ids_per_rank)
    : tiles_{static_cast<unsigned>(config.Mapping.tiles[0]),
             static_cast<unsigned>(config.Mapping.tiles[1]),
             static_cast<unsigned>(config.Mapping.tiles[2])},
      tiles_per_rank_{static_cast<unsigned>(config.Mapping.tilesPerRank[0]),
                      static_cast<unsigned>(config.Mapping.tilesPerRank[1]),
                      static_cast<unsigned>(config.Mapping.tilesPerRank[2])} {
    // Partial-rank samples
    if (tiles_[0] <= tilesPerRank_[0] &&
        tiles_[1] <= tilesPerRank_[1] &&
        tiles_[2] <= tilesPerRank_[2] &&
        tilesPerRank_[0] % tiles_[0] == 0 &&
        tilesPerRank_[1] % tiles_[1] == 0 &&
        tilesPerRank_[2] % tiles_[2] == 0) {
      // The next available rank is partially full ...
      if (next_proc_id > 0) {
        // .. and there's enough space for this sample: place it there.
        // TODO: Handle the case where subsequent samples support different
        // numbers of proc_ids_per_rank.
        if (proc_ids_per_rank() == prev_proc_ids_per_rank &&
            next_proc_id + num_splinters() < prev_proc_ids_per_rank) {
          first_rank_ = next_rank;
          first_proc_id_ = next_proc_id;
        }
        // .. and there isn't enough space: start a new partially-filled rank.
        else {
          first_rank_ = next_rank + 1;
          first_proc_id_ = 0;
        }
      }
      // The next available rank is empty: start a new partially-filled rank.
      else {
        first_rank_ = next_rank;
        first_proc_id_ = 0;
      }
    }
    // Full-rank samples
    else if (tiles_[0] >= tilesPerRank_[0] &&
             tiles_[1] >= tilesPerRank_[1] &&
             tiles_[2] >= tilesPerRank_[2] &&
             tiles_[0] % tilesPerRank_[0] == 0 &&
             tiles_[1] % tilesPerRank_[1] == 0 &&
             tiles_[2] % tilesPerRank_[2] == 0) {
      // If the next available rank is partially full, start from the one
      // after that.
      if (next_proc_id > 0) {
        first_rank_ = next_rank + 1;
      } else {
        first_rank_ = next_rank;
      }
      first_proc_id_ = 0;
    // Invalid samples

    } else {
      std::cerr << "Invalid tiling for sample " << sample_id << std::endl;
      std::cerr.flush();
      exit(1);
    }
  }
  BaseMapping(const BaseMapping& rhs) = delete;
  BaseMapping& operator=(const BaseMapping& rhs) = delete;

public:
  bool isPartialRank() const {
    return
      tiles_[0] <= tiles_per_rank_[0] &&
      tiles_[1] <= tiles_per_rank_[1] &&
      tiles_[2] <= tiles_per_rank_[2];
  }
  ShardID num_shards() const {
    if (isPartialRank()) {
      return 1;
    }
    return
      tiles_[0] / tiles_per_rank_[0] *
      tiles_[1] / tiles_per_rank_[1] *
      tiles_[2] / tiles_per_rank_[2];
  }
  SplinterID num_splinters() const {
    if (isPartialRank()) {
      return tiles_[0] * tiles_[1] * tiles_[2];
    }
    return tiles_per_rank_[0] * tiles_per_rank_[1] * tiles_per_rank_[2];
  }
  AddressSpace next_rank() const {
    if (isPartialRank()) {
      return first_rank_;
    }
    return first_rank_ + num_shards();
  }
  AddressSpace reqd_ranks() const {
    return first_rank_ + num_shards();
  }
  unsigned next_proc_id() const {
    if (isPartialRank()) {
      return first_proc_id_ + num_splinters();
    }
    return 0;
  }
  unsigned proc_ids_per_rank() const {
    return tiles_per_rank_[0] * tiles_per_rank_[1] * tiles_per_rank_[2];
  }
  AddressSpace get_rank(ShardID shard_id) const {
    return first_rank_ + shard_id;
  }
  unsigned get_proc_id(SplinterID splinter_id) const {
    return first_proc_id_ + splinter_id;
  }
  unsigned x_tiles() const {
    return tiles_[0];
  }
  unsigned y_tiles() const {
    return tiles_[1];
  }
  unsigned z_tiles() const {
    return tiles_[2];
  }
  unsigned num_tiles() const {
    return x_tiles() * y_tiles() * z_tiles();
  }

protected:
  unsigned tiles_[3];
  unsigned tiles_per_rank_[3];
  AddressSpace first_rank_;
  unsigned first_proc_id_;
};

#endif // __BASE_MAPPING_H__
