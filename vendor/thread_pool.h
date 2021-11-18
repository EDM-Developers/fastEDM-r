// Adapted from https://github.com/jhasse/ThreadPool/blob/master/ThreadPool.hpp
#pragma once

#include <functional>
#include <future>
#include <queue>
#include <vector>

#ifndef _MSC_VER
#ifndef __APPLE__
#include <pthread.h>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#include <fmt/format.h>
#endif
#endif
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <functional>
#include <future>
#include <queue>
#include <thread>

class ThreadPool
{
public:
  ThreadPool(int);
  template<class F, class... Args>
  auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>;
  ~ThreadPool();
  void set_num_workers(int threads);
  void sync()
  {
    while (tasks.empty() == false) {
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
  }

private:
  // need to keep track of threads so we can join them
  std::vector<std::thread> workers;
  // the task queue
  std::queue<std::function<void()>> tasks;

  // synchronization
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(int threads)
  : stop(false)
{

  set_num_workers(threads);
}

inline void ThreadPool::set_num_workers(int threads)
{
  int numNewThreads = threads - workers.size();

  for (int i = 0; i < numNewThreads; ++i)
    workers.emplace_back([this] {
      for (;;) {
        std::function<void()> task;

        {
          std::unique_lock<std::mutex> lock(this->queue_mutex);
          this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
          if (this->stop && this->tasks.empty())
            return;
          task = std::move(this->tasks.front());
          this->tasks.pop();
        }

        task();
      }
    });
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
{
  using return_type = typename std::result_of<F(Args...)>::type;

  auto task =
    std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

  std::future<return_type> res = task->get_future();
  {
    std::unique_lock<std::mutex> lock(queue_mutex);

    // don't allow enqueueing after stopping the pool
    if (stop)
      throw std::runtime_error("enqueue on stopped ThreadPool");

    tasks.emplace([task]() { (*task)(); });
  }
  condition.notify_one();
  return res;
}

// the destructor joins all threads
ThreadPool::~ThreadPool()
{
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    stop = true;
  }
  condition.notify_all();
  for (std::thread& worker : workers) {
    worker.join();
  }
}

#endif
