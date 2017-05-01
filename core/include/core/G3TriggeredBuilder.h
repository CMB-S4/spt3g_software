#pragma once

#include <G3.h>
#include <G3Frame.h>
#include <G3Logging.h>
#include <G3Module.h>

#include <pthread.h>
#include <ApplePthreadBarrier.h>
#include <deque>
#include <vector>



/**



None of the submodules can use python code since they run in separate threads and the GIL... the fucking GIL


 **/






class G3TriggeredBuilder : public G3Module{
public:
  G3TriggeredBuilder(int is_non_blocking);
  ~G3TriggeredBuilder();
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);


  void AddModule(G3ModulePtr module);  
  void SpawnSubThreads();
  void CleanupThreads();

  //Only 1 thread should call TriggerBlocking
  void Trigger();


private:

  void TriggerBlocking_();
  void TriggerNonBlocking_();
  
  struct ThreadInfoWrapper{
  public:
    ThreadInfoWrapper(G3TriggeredBuilder * context, size_t index) : context(context),
								    index(index){}
    G3TriggeredBuilder * context;
    size_t index;
  };

  std::vector< G3ModulePtr > sub_modules_;
  std::vector< ThreadInfoWrapper > thread_info_;

  std::vector< std::deque< G3FramePtr > > sub_mod_queue_;
  std::vector< G3FramePtr > input_frames_;

  std::deque< G3FramePtr > local_out_queue_;

  std::vector< pthread_t > threads_;
  pthread_barrier_t barrier_wakeup_;
  pthread_barrier_t barrier_done_;

  pthread_mutex_t out_queue_mutex_;

  int should_live_;

  int is_non_blocking_;

  int is_triggering_;
  pthread_mutex_t is_triggering_lock_;
  pthread_barrier_t triggering_barrier_;
  pthread_t trigger_thread_;

  static void *ModuleThreadHelper_(void * thread_info_wrapper);
  static void *TriggerThreadHelper_(void * g3triggered_builder_pntr);


  G3TriggeredBuilder(const G3TriggeredBuilder&); //prevent copy construction      
  G3TriggeredBuilder& operator=(const G3TriggeredBuilder&); //prevent assignment

  SET_LOGGER("G3TriggeredBuilder");
};
G3_POINTERS(G3TriggeredBuilder);
