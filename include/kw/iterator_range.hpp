// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef KW_ITERATOR_RANGE_HPP
#define KW_ITERATOR_RANGE_HPP

#include <boost/range.hpp>

namespace kw {

template <typename Iterator>
boost::iterator_range<Iterator>
make_iterator_range(std::pair<Iterator, Iterator> const& p) {
  return boost::make_iterator_range(p.first, p.second);
}

}  // namespace kw

#endif  // KW_ITERATOR_RANGE_HPP
