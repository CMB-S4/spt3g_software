#include <G3TriggeredBuilder.h>

#include <G3Logging.h>

//log_info
//log_notice
//log_warn
//log_error


G3TriggeredBuilder::~G3TriggeredBuilder(){
  CleanupThreads();
}


G3TriggeredBuilder::G3TriggeredBuilder( int is_non_blocking){
  is_non_blocking_ = is_non_blocking;
  should_live_ = 0;

}

void G3TriggeredBuilder::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){

  out.push_back(frame);
  
  pthread_mutex_lock (&out_queue_mutex_);
  for (auto i= local_out_queue_.begin(); i != local_out_queue_.end(); i++){
    out.push_back(*i);
  }
  pthread_mutex_unlock (&out_queue_mutex_);
}

void G3TriggeredBuilder::AddModule(G3ModulePtr module){
  if (should_live_){
    log_fatal("Added thread when threads are running");
  }
  sub_modules_.push_back(module);
  sub_mod_queue_.push_back(std::deque< G3FramePtr >());
  input_frames_.push_back(G3FramePtr());
}

void G3TriggeredBuilder::Trigger(){
  if (is_non_blocking_) TriggerNonBlocking_();
  else  TriggerBlocking_();
}

void G3TriggeredBuilder::TriggerBlocking_(){
  if (!should_live_){
    log_error("Trying to get values in G3TriggeredBuilder when the child threads are blooming demised");
    return;
  }
  //barrier wakeup
  pthread_barrier_wait( &barrier_wakeup_);
  //barrier done with calc
  pthread_barrier_wait( &barrier_done_);
  
  //collect the queues
  pthread_mutex_lock (&out_queue_mutex_);
  local_out_queue_.clear();
  for (auto i = sub_mod_queue_.begin(); i!= sub_mod_queue_.end(); i++){
    for (auto j = (*i).begin(); j != (*i).end(); j++){
      local_out_queue_.push_back(*j);
    }
  }
  pthread_mutex_unlock (&out_queue_mutex_);
}

void G3TriggeredBuilder::TriggerNonBlocking_(){

  pthread_mutex_lock( &(is_triggering_lock_) );
  if (is_triggering_){
    log_error("Attempting to do non-blocking trigger before the previous trigger has finished running.\n"
	      "Human sacrifice, dogs and cats living together... mass hysteria!");
    pthread_mutex_unlock( &(is_triggering_lock_) );
    return;
  }else {
    is_triggering_ = 1;
    pthread_mutex_unlock( &(is_triggering_lock_) ); 
    pthread_barrier_wait(&(triggering_barrier_));
    return;
  }
}


void G3TriggeredBuilder::SpawnSubThreads(){
  size_t n_threads = sub_modules_.size();
  if (should_live_){
    log_fatal("Attempting to spawn sub threads when threads are running");
  }
  should_live_ = 1;

  pthread_barrier_init( &barrier_wakeup_, NULL, n_threads + 1);
  pthread_barrier_init( &barrier_done_, NULL, n_threads + 1);
  //spawn the pthreads
  threads_ = std::vector< pthread_t > ( n_threads);    
  
  for (size_t i=0; i< n_threads; i++){
    thread_info_.push_back( ThreadInfoWrapper(this, i));
  }

  pthread_mutex_init(&out_queue_mutex_, NULL);
  for (size_t i=0; i< n_threads; i++){
    pthread_create(&(threads_[i]), NULL,
		   &G3TriggeredBuilder::ModuleThreadHelper_,
		   &(thread_info_[i])
		   );
  }

  if (is_non_blocking_){
    is_triggering_ = 0;
    pthread_barrier_init( &triggering_barrier_, NULL, 2);
    pthread_mutex_init(&is_triggering_lock_, NULL);
    pthread_create(&(trigger_thread_), NULL,
		   &G3TriggeredBuilder::TriggerThreadHelper_,
		   this
		   );    
  }

}

void G3TriggeredBuilder::CleanupThreads(){
  if (!should_live_) return;
  should_live_ = 0;
  pthread_barrier_wait( &barrier_wakeup_);
  for (size_t i=0; i < sub_modules_.size(); i++)
    pthread_join(threads_[i], NULL);
  pthread_mutex_destroy(&out_queue_mutex_);

}

void * G3TriggeredBuilder::ModuleThreadHelper_(void * null_tiw){
  G3TriggeredBuilder * tb  = ((ThreadInfoWrapper * )null_tiw)->context;
  size_t index  = ((ThreadInfoWrapper * )null_tiw)->index;
  while (1){
    pthread_barrier_wait( &(tb->barrier_wakeup_));
    if (!tb->should_live_){
      break;
    }
    tb->sub_mod_queue_[index].clear();
    //does thread things
    tb->sub_modules_[index]->Process( tb->input_frames_[index], 
				      tb->sub_mod_queue_[index]  );
    pthread_barrier_wait(&(tb->barrier_done_));
  }
  return NULL;
}


void * G3TriggeredBuilder::TriggerThreadHelper_(void * void_tb){
  G3TriggeredBuilder * tb = (G3TriggeredBuilder *) void_tb;

  // if not triggering trigger

  while (tb->should_live_){
    //have a barrier maybe if not triggering and is living
    pthread_barrier_wait( &(tb->triggering_barrier_));

    if (! tb->should_live_) break;
    tb->TriggerBlocking_();
    
    pthread_mutex_lock( &(tb->is_triggering_lock_) );
    tb->is_triggering_ = 0;
    pthread_mutex_unlock( &(tb->is_triggering_lock_) ); 
  }
  return NULL;
}
