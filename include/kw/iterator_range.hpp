// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef __KW_ITERATOR_RANGE_HPP__
#define __KW_ITERATOR_RANGE_HPP__

#include <boost/range.hpp>

namespace kw {

template <typename Iterator>
boost::iterator_range<Iterator>
make_iterator_range(std::pair<Iterator, Iterator> const& p)
{
    return boost::make_iterator_range(p.first, p.second);
}

}  // namespace kw

#endif  // __KW_ITERATOR_RANGE_HPP__
