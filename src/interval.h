// Copyright 2017 Riku Walve

#ifndef WANDA_INTERVAL_H_
#define WANDA_INTERVAL_H_

class interval_t {
 public:
  size_t left, right;

  interval_t(const size_t _left, const size_t _right) : left(_left), right(_right) {}
  interval_t(const interval_t &interval) : left(interval.left), right(interval.right) {}
  interval_t(interval_t&& interval) noexcept : left(interval.left), right(interval.right) {}

  ~interval_t() noexcept {}

  interval_t& operator=(const interval_t& interval) {
    interval_t tmp(interval);
    *this = std::move(tmp);
    return *this;
  }

  interval_t& operator=(interval_t&& interval) noexcept {
    left = interval.left;
    right = interval.right;
    return *this;
  }

  inline bool operator==(const interval_t &interval) const {
    return (left == interval.left && right == interval.right);
  }

  inline bool operator!=(const interval_t &interval) const {
    return (left != interval.left && right != interval.right);
  }
};

#endif
