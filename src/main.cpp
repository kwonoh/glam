// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <iostream>

#include <boost/geometry/geometries/adapted/boost_array.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/program_options.hpp>
#include <json.hpp>
#include <kw/graph/layout/metric.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_ARRAY_CS(cs::cartesian)
typedef boost::array<double, 2> Point;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    Graph;
typedef
    typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef typename boost::graph_traits<Graph>::edges_size_type edges_size_type;

namespace bpo = boost::program_options;
using json = nlohmann::json;

template <typename PositionMap>
void
load_graph_dot(std::istream& is, Graph& g, PositionMap vp_pos)
{
    auto vp_idx = boost::get(boost::vertex_index, g);
    auto vp_id = boost::make_vector_property_map<std::string>(vp_idx);
    auto vp_x = boost::make_vector_property_map<double>(vp_idx);
    auto vp_y = boost::make_vector_property_map<double>(vp_idx);

    boost::dynamic_properties dp;
    dp.property("vp_id", vp_id);
    dp.property("x", vp_x);
    dp.property("y", vp_y);

    if (!boost::read_graphviz(is, g, dp, "vp_id")) {
        std::cerr << "Unable to load graph" << std::endl;
        return;
    }

    BOOST_FOREACH (auto u, boost::vertices(g)) {
        boost::put(vp_pos, u, Point{boost::get(vp_x, u), boost::get(vp_y, u)});
    }
}

template <typename PositionMap>
void
load_graph_json(std::istream& is, Graph& g, PositionMap vp_pos)
{
    json j;
    is >> j;

    std::vector<vertex_descriptor> vertex_descs;
    for (auto& n : j["nodes"]) {
        auto const u = boost::add_vertex(g);
        vertex_descs.push_back(u);
        boost::put(vp_pos, u, Point{n["x"], n["y"]});
    }

    for (auto& l : j["links"]) {
        std::size_t const s = l["source"];
        std::size_t const t = l["target"];
        boost::add_edge(vertex_descs[s], vertex_descs[t], g);
    }
}

template <typename PositionMap>
void
load_graph(std::string const& filepath, Graph& g, PositionMap vp_pos)
{
    std::cout << "Loading graph: " << filepath << std::endl;

    std::ifstream ifs(filepath);
    if (!ifs.good()) {
        std::cerr << "Unable to open input file: " << filepath << std::endl;
        return;
    }

    std::string const filepath_lc = boost::algorithm::to_lower_copy(filepath);
    if (boost::algorithm::ends_with(filepath_lc, ".dot")) {
        load_graph_dot(ifs, g, vp_pos);
    }
    else if (boost::algorithm::ends_with(filepath_lc, ".json")) {
        load_graph_json(ifs, g, vp_pos);
    }
}

template <typename PositionMap>
void
compute_crosslessness(Graph const& g, PositionMap vp_pos)
{
    edges_size_type num_crossings;
    double const crosslessness = kw::crosslessness(g, vp_pos, num_crossings);
    std::cout << "crosslessness=" << crosslessness
              << " (num_edge_crossings=" << num_crossings << ")" << std::endl;
}

template <typename PositionMap>
void
compute_edge_length_cv(Graph const& g, PositionMap vp_pos)
{
    double cv;
    double norm_cv = kw::edge_length_normalized_cv(g, vp_pos, cv);
    std::cout << "edge_length_cv=" << cv << " (normalized_cv=" << norm_cv << ")"
              << std::endl;
}

template <typename PositionMap>
void
compute_min_angle(Graph const& g, PositionMap vp_pos)
{
    double const value = kw::min_angle_metric(g, vp_pos);
    std::cout << "min_angle=" << value << std::endl;
}

template <typename PositionMap>
void
compute_shape_delaunay(Graph const& g, PositionMap vp_pos)
{
    Graph s;
    auto svp_idx = boost::get(boost::vertex_index, s);
    auto svp_pos = boost::make_vector_property_map<Point>(svp_idx);
    BOOST_FOREACH (auto u, boost::vertices(g)) {
        svp_pos[boost::add_vertex(s)] = vp_pos[u];
    }
    kw::delaunay_triangulation(s, svp_pos);
    double const value = kw::shape_based_metric(g, s);
    std::cout << "shape_delaunay=" << value << std::endl;
}

template <typename PositionMap>
void
compute_shape_gabriel(Graph const& g, PositionMap vp_pos)
{
    Graph s;
    auto svp_idx = boost::get(boost::vertex_index, s);
    auto svp_pos = boost::make_vector_property_map<Point>(svp_idx);
    BOOST_FOREACH (auto u, boost::vertices(g)) {
        svp_pos[boost::add_vertex(s)] = vp_pos[u];
    }
    kw::delaunay_triangulation(s, svp_pos);
    kw::convert_delaunay_triangulation_to_gabriel_graph(s, svp_pos);
    double const value = kw::shape_based_metric(g, s);
    std::cout << "shape_gabriel=" << value << std::endl;
}

template <typename PositionMap>
void
compute_metric(Graph& g, PositionMap vp_pos, std::string const& metric)
{
    std::cout << "Computing metric: " << metric << std::endl;
    if (boost::algorithm::contains(metric, "crosslessness")) {
        compute_crosslessness(g, vp_pos);
    }
    else if (boost::algorithm::contains(metric, "edge_length_cv")) {
        compute_edge_length_cv(g, vp_pos);
    }
    else if (boost::algorithm::contains(metric, "min_angle")) {
        compute_min_angle(g, vp_pos);
    }
    else if (boost::algorithm::contains(metric, "shape_delaunay")) {
        compute_shape_delaunay(g, vp_pos);
    }
    else if (boost::algorithm::contains(metric, "shape_gabriel")) {
        compute_shape_gabriel(g, vp_pos);
    }
}

int
main(int argc, char const* argv[])
{
    bpo::options_description desc("Options");
    desc.add_options()  //
        ("input-file,i", bpo::value<std::vector<std::string>>()->multitoken(),
         "input file(s)")  //
        ("metric,m", bpo::value<std::vector<std::string>>()->multitoken(),
         "metric(s) to compute. Available metrics: "
         "crosslessness, edge_length_cv, shape_gabriel, shape_delaunay")  //
        ("help", "print help message");

    bpo::positional_options_description pod;
    pod.add("input-file", -1);

    bpo::variables_map vm;
    bpo::store(bpo::command_line_parser(argc, argv)
                   .options(desc)
                   .positional(pod)
                   .run(),
               vm);
    bpo::notify(vm);

    if (argc < 1 || vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    if (vm.count("input-file") < 1) {
        std::cerr << "Please provide input file(s)." << std::endl;
        return 1;
    }

    if (vm.count("metric") < 1) {
        std::cerr << "Please provide metric(s) to compute." << std::endl;
        return 1;
    }

    std::vector<std::string> const files =
        vm["input-file"].as<std::vector<std::string>>();
    std::vector<std::string> const metrics =
        vm["metric"].as<std::vector<std::string>>();

    for (std::string const& file : files) {
        Graph g;
        auto vp_idx = boost::get(boost::vertex_index, g);
        auto vp_pos = boost::make_vector_property_map<Point>(vp_idx);
        load_graph(file, g, vp_pos);

        for (std::string const& metric : metrics) {
            compute_metric(g, vp_pos, boost::algorithm::to_lower_copy(metric));
        }
        std::cout << std::endl;
    }

    return 0;
}
