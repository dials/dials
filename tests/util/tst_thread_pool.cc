#include <cassert>
#include <vector>
#include <iostream>
#include <atomic>
#include <dials/util/thread_pool.h>
#include <unistd.h>
#include <cstdlib>

#ifndef _WIN32
#include <sys/wait.h>
#endif

using dials::util::ThreadPool;

void tst_basic_functionality() {
  ThreadPool pool(4);
  std::atomic<int> counter(0);

  // Post some simple tasks
  for (int i = 0; i < 100; ++i) {
    pool.post([&counter]() { counter++; });
  }

  // Wait for all tasks to complete
  pool.wait();

  // Check that all tasks executed
  assert(counter == 100);

  // Test passed
  std::cout << "OK" << std::endl;
}

void tst_concurrent_execution() {
  ThreadPool pool(4);
  std::atomic<int> concurrent_count(0);
  std::atomic<int> max_concurrent(0);
  std::atomic<int> completed(0);

  // Post tasks that simulate work and track concurrency
  for (int i = 0; i < 20; ++i) {
    pool.post([&]() {
      int current = ++concurrent_count;

      // Update max concurrent if needed
      int expected = max_concurrent.load();
      while (expected < current
             && !max_concurrent.compare_exchange_weak(expected, current)) {
        // Keep trying
      }

      // Simulate some work
      for (volatile int j = 0; j < 10000; ++j) {
      }

      --concurrent_count;
      ++completed;
    });
  }

  pool.wait();

  // Should have completed all tasks
  assert(completed == 20);
  // Should have achieved some concurrency (more than 1 thread working)
  assert(max_concurrent > 1);
  // No tasks should be running now
  assert(concurrent_count == 0);

  // Test passed
  std::cout << "OK" << std::endl;
}

void tst_multiple_post_wait_cycles() {
  ThreadPool pool(2);
  std::atomic<int> total_executed(0);

  // Do multiple cycles of post/wait
  for (int cycle = 0; cycle < 5; ++cycle) {
    std::atomic<int> cycle_count(0);

    // Post tasks for this cycle
    for (int i = 0; i < 10; ++i) {
      pool.post([&cycle_count, &total_executed]() {
        cycle_count++;
        total_executed++;
      });
    }

    // Wait for this cycle to complete
    pool.wait();

    // Check this cycle completed
    assert(cycle_count == 10);
  }

  // Check total
  assert(total_executed == 50);

  // Test passed
  std::cout << "OK" << std::endl;
}

void tst_exception_handling() {
  // Just document the current behavior: ThreadPool does not attempt to handle
  // exceptions, instead it just crashes.

#ifdef _WIN32
  std::cout << "Exception handling test skipped on Windows" << std::endl;
  std::cout << "Note: Exceptions thrown in thread pool tasks will terminate the program"
            << std::endl;
  std::cout << "OK" << std::endl;
#else
  // Fork a child process to test exception behavior
  pid_t pid = fork();
  if (pid == 0) {
    // Child process - this should crash
    ThreadPool pool(1);
    pool.post([]() { throw std::runtime_error("Raise an exception"); });
    pool.wait();
    exit(0);  // Should never reach here
  } else {
    // Parent process - wait for child to crash
    int status;
    waitpid(pid, &status, 0);
    assert(!WIFEXITED(status) || WEXITSTATUS(status) != 0);  // Should not exit normally
    std::cout << "OK" << std::endl;
  }
#endif
}

void tst_single_thread_pool() {
  ThreadPool pool(1);
  std::vector<int> execution_order;
  std::atomic<int> counter(0);

  // Post tasks that should execute in order (single thread)
  for (int i = 0; i < 10; ++i) {
    pool.post([&execution_order, &counter, i]() {
      // Simple work simulation
      for (volatile int j = 0; j < 1000; ++j) {
      }
      execution_order.push_back(i);
      counter++;
    });
  }

  pool.wait();

  // All tasks should have executed
  assert(counter == 10);
  assert(execution_order.size() == 10);

  // With single thread, should execute in order
  for (int i = 0; i < 10; ++i) {
    assert(execution_order[i] == i);
  }

  // Test passed
  std::cout << "OK" << std::endl;
}

int main(int argc, char const *argv[]) {
  tst_basic_functionality();
  tst_concurrent_execution();
  tst_multiple_post_wait_cycles();
  tst_exception_handling();
  tst_single_thread_pool();

  return 0;
}
