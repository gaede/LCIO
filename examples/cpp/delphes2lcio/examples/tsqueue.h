/* Copying and distribution of this file, with or without modification,
 * are permitted in any medium without royalty provided the copyright
 * notice and this notice are preserved.  This file is offered as-is,
 * without any warranty. */

#pragma once
#include <utility>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <optional>
#include <cassert>

/* Thread Safe Queue Template, C++17 */
template <typename T>
struct Tsqueue {
/* Create Tsqueue object. Set maximum size of the queue to max_size. */
inline Tsqueue(size_t max_size = -1UL) : maxsize(max_size), end(false) {};

/* Push T to the queue. Many threads can push at the same time.
 * If the queue is full, calling thread will be suspended until
 * some other thread pop() data. */
 void push(const T&);
void push(T&&);

/* Close the queue.
 * Be sure all writing threads done their writes before call this.
  * Push data to closed queue is forbidden. */
void close();

/* Pop and return T from the queue. Many threads can pop at the same time.
 * If the queue is empty, calling thread will be suspended.
  * If the queue is empty and closed, nullopt returned. */
std::optional<T> pop();
private:
std::queue<T> que;
std::mutex mtx;
std::condition_variable cv_empty, cv_full;
const size_t maxsize;
std::atomic<bool> end;
};

/* Usage sample:
#include <thread>
#include <chrono>
#include <iostream>
#include "tsqueue.h"
using namespace std;
Tsqueue<int> q;
void foo()
{
for (int i = 0; i < 4; i++) {
q.push(i);
this_thread::sleep_for(chrono::seconds(1));
}
q.close();
}
int main()
{
thread t(foo);
while(auto x = q.pop())
cout << *x << '\n';
t.join();
}
*/

template<typename T>
void Tsqueue<T>::push(T&& t)
{
std::unique_lock<std::mutex> lck(mtx);
while (que.size() == maxsize && !end)
cv_full.wait(lck);
assert(!end);
que.push(std::move(t));
cv_empty.notify_one();
}

template<typename T>
void Tsqueue<T>::push(const T& t)
{
std::unique_lock<std::mutex> lck(mtx);
while (que.size() == maxsize && !end)
cv_full.wait(lck);
assert(!end);
que.push(std::move(t));
cv_empty.notify_one();
}

template<typename T>
std::optional<T> Tsqueue<T>::pop()
{
std::unique_lock<std::mutex> lck(mtx);
while (que.empty() && !end)
cv_empty.wait(lck);
if (que.empty()) return {};
T t = std::move(que.front());
que.pop();
cv_full.notify_one();
return t;
}

template<typename T>
void Tsqueue<T>::close()
{
end = true;
std::lock_guard<std::mutex> lck(mtx);
cv_empty.notify_all();
cv_full.notify_all();
}

