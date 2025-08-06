#include <cassert>
#include <vector>
#include <iostream>
#include <atomic>
#include <dials/util/thread_pool.h>

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
  ThreadPool pool(2);
  std::atomic<int> successful_tasks(0);
  std::atomic<int> tasks_after_exception(0);

  // Post a mix of normal and exception-throwing tasks
  for (int i = 0; i < 10; ++i) {
    if (i == 5) {
      // Post a task that throws
      pool.post([]() { throw std::runtime_error("Test exception"); });
    } else {
      pool.post([&, i]() {
        successful_tasks++;
        if (i > 5) {
          tasks_after_exception++;
        }
      });
    }
  }

  pool.wait();

  // Should have executed 9 successful tasks
  assert(successful_tasks == 9);
  // Should have executed tasks after the exception
  assert(tasks_after_exception == 4);

  // Test passed
  std::cout << "OK" << std::endl;
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
