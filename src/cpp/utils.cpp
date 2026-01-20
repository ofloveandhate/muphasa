#include "utils.h"



long getMemoryUsage()
{
  struct rusage usage;
  if(0 == getrusage(RUSAGE_SELF, &usage))
    return usage.ru_maxrss; // bytes
  else
    return 0;
}



void mem_usage(double& vm_usage, double& resident_set) {
   vm_usage = 0.0;
   resident_set = 0.0;
   std::ifstream stat_stream("/proc/self/stat",std::ios_base::in); //get info from proc directory
   //create some variables to get info
   std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
   stat_stream.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured to use 2MB pages
   vm_usage = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

void get_mem_usage(double& vm_usage, double& resident_set){

  /* // Commented for Linuxed -- MB
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS == task_info(mach_task_self(),
                                  TASK_BASIC_INFO, (task_info_t)&t_info,
                                  &t_info_count))
    {
        printf("Proccess :: virtual memory: %lu Mb\n resident memory: %lu Mb \n", (t_info.virtual_size/1024)/1024, (t_info.resident_size/1024)/1024);
        vm_usage = ((double)t_info.virtual_size/1024)/1024;
        resident_set = ((double)t_info.resident_size/1024)/1024;
    }
  */
}

