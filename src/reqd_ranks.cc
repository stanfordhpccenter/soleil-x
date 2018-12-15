#include <iostream>

#include "config_schema.h"
#include "base_mapping.h"

#define EQUALS(s1, s2) (strcmp((s1), (s2)) == 0)

int main(int argc, char* argv[]) {
  unsigned num_samples = 0;
  unsigned reqd_ranks = 0;
  AddressSpace next_rank = 0;
  unsigned next_proc_id = 0;
  unsigned prev_proc_ids_per_rank = 1;
  auto process_config = [&](const Config& config) {
    BaseMapping mapping(config, num_samples,
                        next_rank, next_proc_id,
                        prev_proc_ids_per_rank);
    num_samples++;
    reqd_ranks = mapping.reqd_ranks();
    next_rank = mapping.next_rank();
    next_proc_id = mapping.next_proc_id();
    prev_proc_ids_per_rank = mapping.proc_ids_per_rank();
  };
  for (int i = 0; i < argc; ++i) {
    if (EQUALS(argv[i], "-i") && i < argc-1) {
      Config config;
      parse_Config(&config, argv[i+1]);
      process_config(config);
    } else if (EQUALS(argv[i], "-m") && i < argc-1) {
      MultiConfig mc;
      parse_MultiConfig(&mc, argv[i+1]);
      process_config(mc.configs[0]);
      process_config(mc.configs[1]);
    }
  }
  std::cout << reqd_ranks << std::endl;
  return 0;
}
