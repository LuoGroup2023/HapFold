#ifndef _BUBBLE_CHAIN_H_
#define _BUBBLE_CHAIN_H_
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <set>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "graph.h"
#include <map>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <queue>
#include <numeric>
#include <cmath>
#include <stack>
#include <functional>
using namespace std;

set<uint32_t> complex_components;
typedef struct
{
    uint32_t p;             // the optimal parent vertex
    uint32_t d;             // the shortest distance from the initial vertex
    uint32_t c;             // max count of positive reads
    uint32_t m;             // max count of negative reads
    uint32_t np;            // max count of non-positive reads
    uint32_t nc;            // max count of reads, no matter positive or negative
    uint32_t r : 31, s : 1; // r: the number of remaining incoming arc; s: state
                            // s: state, s=0, this edge has not been visited, otherwise, s=1
} binfo_t;

typedef struct
{
    /// all information for each node
    binfo_t *a;
    kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
    kvec_t(uint32_t) T; // set of tips
    kvec_t(uint32_t) b; // visited vertices
    kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;

typedef struct
{
    uint32_t id;
    uint32_t begNode;
    uint32_t endNode;
    vector<vector<asg_arc_t *>> valid_paths;
    vector<vector<asg_arc_t *>> paths;
    vector<set<uint32_t>> paths_nodes;
    set<uint32_t> starting_arcs;
    set<uint32_t> ending_arcs;
} bubble_t;

typedef struct bubble_chain_t
{
    uint32_t id;
    vector<bubble_t *> bubbles; // bubbles in order
    int *usl;                   // unique sequence length of unique sequences between bubbles
};

typedef struct bubble_chain_graph_t
{
    vector<uint32_t> normal_nodes;
    vector<bubble_chain_t *> condensed_bubble_chain;
};

#define NORMAL_NODE 0;
#define BUBBLE_END_BEGIN 1;
#define BUBBLE_INSIDE 2;
/*
 * Get the sources of a single component
 */
vector<uint32_t> get_sources(asg_t *g)
{
    cout << "start source search" << endl;
    vector<uint32_t> sources;

    uint32_t n_vtx = g->n_seq * 2;
    bool frontier_added[n_vtx];
    for (int i = 0; i < n_vtx; i++)
    {
        frontier_added[i] = false;
    }
    vector<uint32_t> frontier;
    uint32_t s0 = 0; // arbitrary vertex
    frontier.push_back(s0);
    frontier_added[s0] = true;
    while (!frontier.empty())
    {
        uint32_t u = frontier.back();
        frontier.pop_back();
        // cout << u << "\t" << g->seq[u/2].name << endl;

        uint32_t num_incoming_arcs = asg_arc_n(g, u ^ 1); // 1
        if (num_incoming_arcs == 0)
        { // check if curr_v is a source by checking if no parents
            sources.push_back(u);
        }

        // add parents to frontier
        // uint32_t num_incoming_arcs = asg_arc_n(g, u^1);  // 1
        asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u ^ 1); // incoming_arcs_complement[0].v^1 = 22;
        for (int vi = 0; vi < num_incoming_arcs; vi++)
        {
            int p = incoming_arcs_complement[vi].v ^ 1;
            if (!frontier_added[p])
            {
                frontier.push_back(p);
                frontier_added[p] = true;
            }
        }

        // add children to frontier
        uint32_t num_outgoing_arcs = asg_arc_n(g, u); // 2
        asg_arc_t *outgoing_arcs = asg_arc_a(g, u);   // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
        for (int vi = 0; vi < num_outgoing_arcs; vi++)
        {
            int c = outgoing_arcs[vi].v;
            if (!frontier_added[c])
            {
                frontier.push_back(c);
                frontier_added[c] = true;
            }
        }
    }
    cout << "finish source search" << endl;
    return sources;
}

vector<uint32_t> get_topological_order(asg_t *g)
{
    cout << "start get topological order" << endl;

    vector<uint32_t> reverse_topological_order;

    uint32_t n_vtx = g->n_seq * 2;
    bool frontier_added[n_vtx], child_frontier_added[n_vtx];
    for (int u = 0; u < n_vtx; u++)
    {
        frontier_added[u] = false;
        child_frontier_added[u] = false;
    }
    vector<uint32_t> component, frontier, child_frontier, child_vi_frontier;
    uint32_t s0 = 0; // arbitrary vertex
    frontier.push_back(s0);
    frontier_added[s0] = true;
    while (!frontier.empty())
    {
        // for (int u=0; u<frontier.size(); u++) {
        //     cout << frontier[u] << ", ";
        // }
        // cout << endl;

        uint32_t u = frontier.back();
        frontier.pop_back();
        if (!child_frontier_added[u])
        {
            child_frontier.push_back(u);
            child_vi_frontier.push_back(0);
            child_frontier_added[u] = true;
        }
        while (!child_frontier.empty())
        {
            // for (int u=0; u<child_frontier.size(); u++) {
            //     cout << child_frontier[u] << ", ";
            // }
            // cout << endl;

            uint32_t u = child_frontier.back();
            uint32_t vi = child_vi_frontier.back();

            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            if (vi < num_outgoing_arcs)
            {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
                uint32_t v = outgoing_arcs[vi].v;
                child_vi_frontier.back()++;

                if (!child_frontier_added[v])
                {
                    child_frontier.push_back(v);
                    child_vi_frontier.push_back(0);
                    frontier_added[v] = true;
                    child_frontier_added[v] = true;
                }
            }
            else
            {
                int num_incoming_arcs = asg_arc_n(g, u ^ 1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u ^ 1);
                for (int uivj = 0; uivj < num_incoming_arcs; uivj++)
                {
                    int vj = incoming_arcs_complement[uivj].v ^ 1;
                    if (!frontier_added[vj])
                    {
                        frontier.push_back(vj);
                        frontier_added[vj] = true;
                    }
                }
                child_frontier.pop_back();
                child_vi_frontier.pop_back();
                reverse_topological_order.push_back(u);
            }
        }
    }

    reverse(reverse_topological_order.begin(), reverse_topological_order.end());
    vector<uint32_t> topological_order = reverse_topological_order; // TODO
    // for (int u=0; u<topological_order.size(); u++) {
    //     cout << u << "\t" << g->seq[topological_order[u]/2].name << endl;
    // }

    cout << "finish get topological order" << endl;
    return topological_order;
}

vector<uint32_t> get_bubble_ends(asg_t *g, vector<uint32_t> sources, vector<uint32_t> *bubble_beginnings, vector<uint32_t> *bubble_ends)
{
    cout << "start get bubble ends" << endl;

    uint32_t n_vtx = g->n_seq * 2;
    bool seen[n_vtx], visited[n_vtx];
    for (int u = 0; u < n_vtx; u++)
    {
        seen[u] = false;
        visited[u] = false;
    }
    for (int ui = 0; ui < sources.size(); ui++)
    {
        cout << "Source: " << sources[ui] << "\t" << g->seq[sources[ui] / 2].name << endl;
    }
    vector<uint32_t> frontier = sources;
    int num_seen = frontier.size(); // seen = seen but not visited
    while (!frontier.empty())
    {
        uint32_t u = frontier.back();
        frontier.pop_back();
        // cout << u << "\t" << g->seq[u/2].name << endl;  // topological order
        if (not visited[u])
        {
            visited[u] = true;
            seen[u] = false;
            num_seen--;

            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            uint32_t num_incoming_arcs = asg_arc_n(g, u ^ 1);
            if (num_seen == 0)
            {
                if (num_outgoing_arcs > 1)
                {
                    bubble_beginnings->push_back(u); // last one is meaningless/unpaired
                    cout << "Beginning:" << u << "\t" << g->seq[u / 2].name << endl;
                }
                bool sources_contains_u = find(sources.begin(), sources.end(), u) != sources.end();
                if (num_incoming_arcs > 1)
                {
                    bubble_ends->push_back(u); // first one is meaningless/unpaired
                    cout << "Ending:" << u << "\t" << g->seq[u / 2].name << endl;
                }
            }

            // add children to frontier if they have no other unvisited parents
            // uint32_t num_outgoing_arcs = asg_arc_n(g, u);  // 2
            asg_arc_t *outgoing_arcs = asg_arc_a(g, u); // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
            for (int vi = 0; vi < num_outgoing_arcs; vi++)
            {
                int c = outgoing_arcs[vi].v;

                uint32_t num_unvisited_wives = 0;
                // uint32_t num_incoming_arcs = asg_arc_n(g, c^1);
                num_incoming_arcs = asg_arc_n(g, c ^ 1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, c ^ 1);
                for (int viuj = 0; viuj < num_incoming_arcs; viuj++)
                {
                    int w = incoming_arcs_complement[viuj].v ^ 1; // wives
                    if (!visited[w])
                    {
                        num_unvisited_wives++;
                    }
                }

                if (!seen[c])
                {
                    num_seen++;
                }
                seen[c] = true;
                if (num_unvisited_wives == 0)
                {
                    frontier.push_back(c);
                }
            }
        }
    }
    cout << "finish get bubble ends: " << num_seen << endl;
    for (int u = 0; u < n_vtx; u++)
    {
        if (!visited[u])
        {
            cout << "Unvisited: " << u << endl;
        }
    }
    return *bubble_ends;
}

bubble_t *detect_bubble(asg_t *g, uint32_t source)
{
    uint32_t num_source_outgoing_arcs = asg_arc_n(g, source);
    uint32_t n_vtx = g->n_seq * 2;
    uint32_t num_seen = 1;
    if (num_source_outgoing_arcs < 2)
    {
        return nullptr;
    }

    bool seen[n_vtx], visited[n_vtx];
    for (int u = 0; u < n_vtx; u++)
    {
        seen[u] = false;
        visited[u] = false;
    }
    vector<uint32_t> frontier;

    frontier.push_back(source);
    seen[source] = true;

    while (!frontier.empty())
    {
        uint32_t v = frontier.back();
        frontier.pop_back();
        // cout << u << "\t" << g->seq[v/2].name << endl;  // topological order
        assert(seen[v]);
        assert(!visited[v]);
        visited[v] = true;
        seen[v] = false;
        num_seen--;
        uint32_t num_outgoing_arcs = asg_arc_n(g, v);

        if (num_outgoing_arcs == 0)
        {
            return nullptr;
        }

        // add children to frontier if they have no other unvisited parents
        // uint32_t num_outgoing_arcs = asg_arc_n(g, v);  // 2
        asg_arc_t *outgoing_arcs = asg_arc_a(g, v); // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
        set<uint32_t> edge_set;
        for (int vi = 0; vi < num_outgoing_arcs; vi++)
        {
            uint32_t u = outgoing_arcs[vi].v;
            if (edge_set.insert(u).second)
            {
                if (u == v || u == (v ^ 1) || visited[u ^ 1] || u == source)
                {
                    return nullptr;
                }

                assert(!visited[u]);
                if (!seen[u])
                {
                    num_seen++;
                }
                seen[u] = true;
                assert(asg_arc_n(g, u ^ 1) >= 1);
                bool has_unvisited_parents = false;
                uint32_t num_incoming_arcs = asg_arc_n(g, u ^ 1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u ^ 1);
                for (int viuj = 0; viuj < num_incoming_arcs; viuj++)
                {
                    uint32_t p = (incoming_arcs_complement[viuj].v) ^ 1; // parent
                    if (!visited[p])
                        has_unvisited_parents = true;
                }

                if (!has_unvisited_parents)
                {
                    frontier.push_back(u);
                }
            }
        }

        if (frontier.size() == 1 && num_seen == 1 && seen[frontier.front()])
        {
            uint32_t t = frontier.back();
            frontier.pop_back();
            uint32_t num_outgoing_arcs = asg_arc_n(g, t);
            if (num_outgoing_arcs > 0)
            {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, t);
                for (int vi = 0; vi < num_outgoing_arcs; vi++)
                {
                    if (outgoing_arcs[vi].v == source)
                    {
                        return nullptr;
                    }
                }
            }
            bubble_t *result = new bubble_t();
            result->begNode = source;
            result->endNode = t;
            return result;
        }
    }
    return nullptr;
}

map<uint32_t, map<uint32_t, set<uint32_t>>> *get_bubble_chain_graph(asg_t *g, set<uint32_t> bubble_chain_begin_end, string output_directory)
{
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes = new map<uint32_t, map<uint32_t, set<uint32_t>>>();
    for (auto bubble_chain_begin : bubble_chain_begin_end)
    {
        map<uint32_t, set<uint32_t>> current_map;
        vector<asg_arc_t *> arc_stack;
        vector<uint32_t> node_stack; // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_chain_begin);
        node_vi_stack.push_back(0);
        bool visited[g->n_seq * 2];
        for (int u = 0; u < g->n_seq * 2; u++)
        {
            visited[u] = false;
        }
        while (!node_stack.empty())
        {
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
            visited[u] = true;
            int flag = 0;

            if (u != bubble_chain_begin && find(bubble_chain_begin_end.begin(), bubble_chain_begin_end.end(), u) != bubble_chain_begin_end.end())
            {
                flag = 1;
            }
            else
            {
                for (auto b : current_map)
                {
                    if (flag)
                    {
                        break;
                    }
                    for (int i = 0; i < num_outgoing_arcs; i++)
                    {
                        if (b.second.find(outgoing_arcs[i].v) != b.second.end())
                        {
                            u = b.first;
                            flag = 2;
                            break;
                        }
                    }
                }
            }

            if (flag == 1)
            {
                if (current_map.find(u) == current_map.end())
                {
                    set<uint32_t> new_set;
                    current_map.insert({u, new_set});
                }
                if (arc_stack.size() > 0)
                {
                    for (int a = 0; a < arc_stack.size() - 1; a++)
                    {
                        current_map[u].insert(arc_stack[a]->v);
                    }
                }
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }
            else if (flag == 2)
            {
                if (current_map.find(u) == current_map.end())
                {
                    set<uint32_t> new_set;
                    current_map.insert({u, new_set});
                }
                if (arc_stack.size() > 0)
                {
                    for (int a = 0; a < arc_stack.size(); a++)
                    {
                        current_map[u].insert(arc_stack[a]->v);
                    }
                }
            }

            if (vi < num_outgoing_arcs)
            {
                while (vi < num_outgoing_arcs)
                {
                    uint32_t v = outgoing_arcs[vi].v;
                    node_vi_stack.back()++;
                    if (!visited[v])
                    {
                        arc_stack.push_back(outgoing_arcs + vi);
                        node_stack.push_back(v);
                        node_vi_stack.push_back(0);
                        break;
                    }
                    else
                    {
                        vi++;
                    }
                }
            }
            if (vi >= num_outgoing_arcs)
            {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
        }
        bubble_chain_begin_end_nodes->insert({bubble_chain_begin, current_map});
    }

    for (auto a : *bubble_chain_begin_end_nodes)
    {
        for (auto b : a.second)
        {
            (*bubble_chain_begin_end_nodes)[a.first][b.first].insert(a.first);
            (*bubble_chain_begin_end_nodes)[a.first][b.first].insert(b.first);
        }
    }
    system((string("rm -r ") + output_directory).c_str());
    mkdir(output_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    ofstream main_graph(output_directory + string("/") + string("main_graph.gfa"));
    for (auto a : *bubble_chain_begin_end_nodes)
    {
        if (a.first % 2 == 0)
        {
            main_graph << "S\t" << g->seq[a.first >> 1].name << "\t" << "*" << endl;
        }
    }
    for (auto a : *bubble_chain_begin_end_nodes)
    {
        for (auto b : a.second)
        {
            main_graph << "L\t" << g->seq[a.first >> 1].name << "\t" << (a.first % 2 == 0 ? "+" : "-") << "\t" << g->seq[b.first >> 1].name << "\t" << (b.first % 2 == 0 ? "+" : "-") << "\t0M\t" << endl;
        }
    }
    main_graph.close();

    ofstream non_length_graph(output_directory + string("/") + string("no_longth_full_graph.gfa"));
    for (int i = 0; i < g->n_seq; i++)
    {
        non_length_graph << "S\t" << g->seq[i].name << "\t" << "*" << endl;
    }

    for (int i = 0; i < g->n_seq * 2; i++)
    {
        uint32_t num_outgoing_arcs = asg_arc_n(g, i);
        asg_arc_t *outgoing_arcs = asg_arc_a(g, i);
        for (int q = 0; q < num_outgoing_arcs; q++)
        {
            non_length_graph << "L\t" << g->seq[i >> 1].name << "\t" << (i % 2 == 0 ? "+" : "-") << "\t" << g->seq[outgoing_arcs[q].v >> 1].name << "\t" << (outgoing_arcs[q].v % 2 == 0 ? "+" : "-") << "\t0M\t" << endl;
        }
    }
    non_length_graph.close();

    cout << "finish get bubble chain" << endl;

    return bubble_chain_begin_end_nodes;
}

map<uint32_t, map<uint32_t, set<uint32_t>>> *get_bubbles(asg_t *g, string output_directory, uint32_t **connection_forward, uint32_t **connection_backward, int detect_cycle_flag = false, int detect_branch_flag = false)
{
    uint32_t **connections_count;
    CALLOC(connections_count, g->n_seq);
    for (int i = 0; i < g->n_seq; i++)
    {
        CALLOC(connections_count[i], g->n_seq);
        memset(connections_count[i], 0, sizeof(*connections_count[i]));
    }
    for (int i = 0; i < g->n_seq; i++)
    {
        for (int j = 0; j < g->n_seq; j++)
        {
            connections_count[i][j] = connection_forward[i][j] + connection_backward[i][j];
            connections_count[i][j] += connection_forward[j][i] + connection_backward[j][i];
        }
    }
    cout << "start get bubbles" << endl;

    uint32_t n_vtx = g->n_seq * 2;
    int node_type[n_vtx];

    vector<bubble_t *> bubble_by_ending_begining[n_vtx];
    vector<bubble_t *> bubbles;
    for (uint32_t i = 0; i < n_vtx; i++)
    {
        node_type[i] = 0;
        bubble_t *result = detect_bubble(g, i);
        if (result != nullptr)
        {
            result->id = i;
            bubble_by_ending_begining[result->begNode].push_back(result);
            bubble_by_ending_begining[result->endNode].push_back(result);
            bubbles.push_back(result);
        }
    }

    for (int b = 0; b < bubbles.size(); b++)
    {
        bubble_t *bubble = bubbles[b];
        uint32_t bubble_beginning = bubble->begNode;
        uint32_t bubble_end = bubble->endNode;

        // cout << "start get bubble paths from " << g->seq[bubble_beginning /2].name<< " to " << g->seq[bubble_end /2].name << endl;

        vector<asg_arc_t *> arc_stack;
        vector<uint32_t> node_stack; // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_beginning);
        node_vi_stack.push_back(0);
        int stack_count = 0;
        uint32_t lastNode = bubble_beginning;
        while (!node_stack.empty())
        {
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            if (u == bubble_end)
            {
                for (int ui = 0; ui < node_stack.size(); ui++)
                {
                    // cout << g->seq[node_stack[ui]/2].name << " ";
                    if (ui == 0 || ui == node_stack.size() - 1)
                    {
                        if (node_type[node_stack[ui]] != 2)
                        {
                            node_type[node_stack[ui]] = BUBBLE_END_BEGIN;
                        }
                    }
                    else
                    {
                        node_type[node_stack[ui]] = BUBBLE_INSIDE;
                    }
                }
                stack_count++;
                // cout << endl;

                bubble->paths.push_back(arc_stack);
                bubble->starting_arcs.insert(arc_stack.front()->v);
                bubble->ending_arcs.insert(lastNode);
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }

            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            if (vi < num_outgoing_arcs)
            {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
                uint32_t v = outgoing_arcs[vi].v;
                node_vi_stack.back()++;
                arc_stack.push_back(outgoing_arcs + vi);
                node_stack.push_back(v);
                node_vi_stack.push_back(0);
            }
            else
            {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
            lastNode = u;
        }
        // cout << stack_count << endl;
        // cout << "finish get bubble paths from " << bubble_beginning << " to " << bubble_end << endl;
    }

    set<uint32_t> bubble_chain_end_begin;
    map<uint32_t, uint32_t> pure_outgoing_num;
    map<uint32_t, uint32_t> pure_incoming_num;
    set<uint32_t> out_only;
    for (uint32_t i = 0; i < n_vtx; i++)
    {
        if (node_type[i] != 2)
        {
            int num_outgoing_arcs = asg_arc_n(g, i);
            int num_incoming_arcs = asg_arc_n(g, i ^ 1);
            asg_arc_t *all_outgoing_arcs = asg_arc_a(g, i);
            asg_arc_t *all_incoming_arcs = asg_arc_a(g, i ^ 1);
            set<uint32_t> valid_outgoing_arcs, valid_incoming_arcs;
            for (int c = 0; c < num_outgoing_arcs; c++)
            {
                valid_outgoing_arcs.insert(all_outgoing_arcs[c].v);
            }
            for (int c = 0; c < num_incoming_arcs; c++)
            {
                valid_incoming_arcs.insert(all_incoming_arcs[c].v ^ 1);
            }
            num_incoming_arcs = valid_incoming_arcs.size();
            num_outgoing_arcs = valid_outgoing_arcs.size();
            if (node_type[i] == 1)
            {
                vector<bubble_t *> bubbles = bubble_by_ending_begining[i];
                set<uint32_t> outgoing_arcs_set, incoming_arcs_set;
                sort(bubbles.begin(), bubbles.end(), [](const auto &lhs, const auto &rhs)
                     { return lhs->starting_arcs.size() > rhs->starting_arcs.size(); });
                for (int q = 0; q < bubbles.size(); q++)
                {
                    bubble_t *bubble = bubbles[q];
                    if (bubble->begNode == i)
                    {
                        bool flag = false;
                        // cout << g->seq[bubble->begNode/2].name << " to " << g->seq[bubble->endNode/2].name << " starts: " <<bubble->starting_arcs.size() << " ends: " << bubble->ending_arcs.size() <<endl;
                        for (auto c : bubble->starting_arcs)
                        {
                            if (outgoing_arcs_set.insert(c).second)
                            {
                                num_outgoing_arcs--;
                                flag = true;
                            }
                        }
                        if (flag)
                            num_outgoing_arcs++;
                    }
                }
                sort(bubbles.begin(), bubbles.end(), [](const auto &lhs, const auto &rhs)
                     { return lhs->ending_arcs.size() > rhs->ending_arcs.size(); });
                for (int q = 0; q < bubbles.size(); q++)
                {
                    bubble_t *bubble = bubbles[q];
                    if (bubble->endNode == i)
                    {
                        bool flag = false;
                        for (auto c : bubble->ending_arcs)
                        {
                            if (incoming_arcs_set.insert(c).second)
                            {
                                num_incoming_arcs--;
                                flag = true;
                            }
                        }
                        if (flag)
                            num_incoming_arcs++;
                    }
                }
            }
            pure_outgoing_num[i] = num_outgoing_arcs;
            pure_incoming_num[i] = num_incoming_arcs;
            if (num_incoming_arcs == 0 && num_outgoing_arcs > 0)
            {
                out_only.insert(i);
            }
            if (!(num_incoming_arcs == 1 && num_outgoing_arcs == 1))
            {
                // cout << g->seq[i/2].name << ": incoming " << num_incoming_arcs << "; outgoing " << num_outgoing_arcs <<"; bubbles:" << bubble_by_ending_begining[i].size() <<endl;
                bubble_chain_end_begin.insert(i);
            }
            // TODO:: implement
        }
    }

    // set<uint32_t> not_accessible;
    // for(auto beg: out_only){
    //     vector<asg_arc_t*> arc_stack;
    //     vector<uint32_t> node_stack;  // DFS
    //     vector<uint32_t> node_vi_stack;
    //     node_stack.push_back(beg);
    //     node_vi_stack.push_back(0);
    //     bool visited[g->n_seq * 2];
    //     for (int u=0; u<g->n_seq * 2; u++) {
    //         visited[u] = false;
    //     }
    //     uint32_t cur_node = beg;
    //     uint32_t linear_len = 1;
    //     while(pure_outgoing_num.find(cur_node) != pure_outgoing_num.end() && pure_outgoing_num[cur_node]!=1) {
    //         int num_outgoing_arcs = asg_arc_n(g, i);
    //         asg_arc_t* all_outgoing_arcs = asg_arc_a(g, i);

    //     }
    // }

    // set<string> names;
    // for(auto c : bubble_chain_end_begin){
    //     // cout << c << endl;
    //     names.insert(g->seq[c/2].name);
    // }

    // for(auto c : names){
    //     cout << c << endl;
    // }
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes_buf = get_bubble_chain_graph(g, bubble_chain_end_begin, output_directory);
    for (auto a : *bubble_chain_begin_end_nodes_buf)
    {
        for (auto b : a.second)
        {
            for (auto bub : bubble_by_ending_begining[a.first])
            {
                if (bub->begNode == a.first)
                {
                    for (auto node : bub->starting_arcs)
                    {
                        (*bubble_chain_begin_end_nodes_buf)[a.first][b.first].insert(node);
                    }
                }
            }
            for (auto bub : bubble_by_ending_begining[b.first])
            {
                if (bub->endNode == b.first)
                {
                    for (auto node : bub->ending_arcs)
                    {
                        (*bubble_chain_begin_end_nodes_buf)[a.first][b.first].insert(node);
                    }
                }
            }
        }
    }
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes_buf2 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();
    for (auto a : *bubble_chain_begin_end_nodes_buf)
    {
        for (auto b : a.second)
        {
            if ((*bubble_chain_begin_end_nodes_buf2).find(a.first) == (*bubble_chain_begin_end_nodes_buf2).end())
            {
                (*bubble_chain_begin_end_nodes_buf2)[a.first] = map<uint32_t, set<uint32_t>>();
            }
            if ((*bubble_chain_begin_end_nodes_buf2).find(b.first ^ 1) == (*bubble_chain_begin_end_nodes_buf2).end())
            {
                (*bubble_chain_begin_end_nodes_buf2)[b.first ^ 1] = map<uint32_t, set<uint32_t>>();
            }
            (*bubble_chain_begin_end_nodes_buf2)[a.first][b.first] = b.second;
            (*bubble_chain_begin_end_nodes_buf2)[b.first ^ 1][a.first ^ 1] = b.second;
        }
    }
    delete bubble_chain_begin_end_nodes_buf;
    set<uint32_t> unaccessable;
    set<uint32_t> not_filtered;
    for (auto a : *bubble_chain_begin_end_nodes_buf2)
    {
        if (pure_outgoing_num[a.first] > 1 && pure_incoming_num[a.first] > 0)
        {
            map<uint32_t, set<uint32_t>> non_outgoing;
            for (auto b : a.second)
            {
                if (pure_outgoing_num[b.first] == 0)
                {
                    non_outgoing[b.first] = b.second;
                }
            }
            if (non_outgoing.size() == 1)
            {
                for (auto x : non_outgoing)
                {
                    if (x.second.size() <= 8)
                    {
                        unaccessable.insert(x.first >> 1);
                    }
                    else
                    {
                        not_filtered.insert(x.first >> 1);
                    }
                }
            }
            else
            {
                uint32_t greater_count = 0;
                for (auto x : non_outgoing)
                {
                    if (x.second.size() > 8)
                    {
                        greater_count++;
                    }
                }
                if (greater_count > 0)
                {
                    for (auto x : non_outgoing)
                    {
                        if (x.second.size() <= 8)
                        {
                            unaccessable.insert(x.first >> 1);
                        }
                        else
                        {
                            not_filtered.insert(x.first >> 1);
                        }
                    }
                }
                else
                {
                    uint32_t max_idx = -1;
                    uint32_t max_count = 0;
                    for (auto x : non_outgoing)
                    {
                        uint32_t cur_len = 0;
                        for (auto i : x.second)
                        {
                            cur_len += g->seq[i >> 1].len;
                        }
                        if (cur_len >= max_count)
                        {
                            max_count = cur_len;
                            max_idx = x.first;
                        }
                    }
                    for (auto x : non_outgoing)
                    {
                        if (x.first != max_idx)
                        {
                            unaccessable.insert(x.first >> 1);
                        }
                        else
                        {
                            not_filtered.insert(x.first >> 1);
                        }
                    }
                }
            }
        }
    }
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes_buf3 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_end_begin_nodes_buf3 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();

    for (auto a : *bubble_chain_begin_end_nodes_buf2)
    {
        if (unaccessable.find(a.first >> 1) == unaccessable.end())
        {
            for (auto b : a.second)
            {
                if (unaccessable.find(b.first >> 1) == unaccessable.end() && b.second.size() > 1)
                {
                    if ((*bubble_chain_begin_end_nodes_buf3).find(a.first) == (*bubble_chain_begin_end_nodes_buf3).end())
                    {
                        (*bubble_chain_begin_end_nodes_buf3)[a.first] = map<uint32_t, set<uint32_t>>();
                    }
                    if ((*bubble_chain_begin_end_nodes_buf3).find(b.first ^ 1) == (*bubble_chain_begin_end_nodes_buf3).end())
                    {
                        (*bubble_chain_begin_end_nodes_buf3)[b.first ^ 1] = map<uint32_t, set<uint32_t>>();
                    }
                    (*bubble_chain_begin_end_nodes_buf3)[a.first][b.first] = b.second;
                    (*bubble_chain_begin_end_nodes_buf3)[b.first ^ 1][a.first ^ 1] = b.second;
                    if ((*bubble_chain_end_begin_nodes_buf3).find(b.first) == (*bubble_chain_end_begin_nodes_buf3).end())
                    {
                        (*bubble_chain_end_begin_nodes_buf3)[b.first] = map<uint32_t, set<uint32_t>>();
                    }
                    if ((*bubble_chain_end_begin_nodes_buf3).find(a.first ^ 1) == (*bubble_chain_end_begin_nodes_buf3).end())
                    {
                        (*bubble_chain_end_begin_nodes_buf3)[a.first ^ 1] = map<uint32_t, set<uint32_t>>();
                    }
                    (*bubble_chain_end_begin_nodes_buf3)[b.first][a.first] = b.second;
                    (*bubble_chain_end_begin_nodes_buf3)[a.first ^ 1][b.first ^ 1] = b.second;
                }
            }
        }
    }
    delete bubble_chain_begin_end_nodes_buf2;
    set<vector<uint32_t>> extensions;
    cout << "Start connect single branches" << endl;
    for (auto beg : (*bubble_chain_begin_end_nodes_buf3))
    {
        for (auto end : beg.second)
        {
            vector<uint32_t> to_expand;
            set<uint32_t> visited;
            visited.insert(beg.first);
            visited.insert(end.first);
            to_expand.push_back(beg.first);
            to_expand.push_back(end.first);
            while ((*bubble_chain_begin_end_nodes_buf3).find(to_expand[to_expand.size() - 1]) != (*bubble_chain_begin_end_nodes_buf3).end() && (*bubble_chain_begin_end_nodes_buf3)[to_expand[to_expand.size() - 1]].size() == 1 && (*bubble_chain_end_begin_nodes_buf3).find(to_expand[to_expand.size() - 1]) != (*bubble_chain_end_begin_nodes_buf3).end() && (*bubble_chain_end_begin_nodes_buf3)[to_expand[to_expand.size() - 1]].size() == 1)
            {
                bool seen = false;
                ;
                for (auto i : (*bubble_chain_begin_end_nodes_buf3)[to_expand[to_expand.size() - 1]])
                {
                    if (visited.find(i.first) != visited.end())
                    {
                        seen = true;
                        break;
                    }
                    // cout << i.first << endl;
                    to_expand.push_back(i.first);
                    visited.insert(i.first);
                }
                if (seen)
                {
                    break;
                }
            }
            extensions.insert(to_expand);
        }
    }
    cout << "Start filter short connections:\t" << extensions.size() << endl;
    vector<vector<uint32_t>> pathes;
    for (auto path : extensions)
    {
        set<uint32_t> path_nodes;
        for (auto n : path)
        {
            path_nodes.insert(n);
        }
        bool valid = true;
        for (auto to_check : extensions)
        {
            set<uint32_t> to_check_nodes;
            for (auto n : to_check)
            {
                to_check_nodes.insert(n);
            }
            bool found = true;
            for (auto n : path_nodes)
            {
                if (to_check_nodes.find(n) == to_check_nodes.end())
                {
                    found = false;
                    break;
                }
            }
            if (found && to_check_nodes.size() > path_nodes.size())
            {
                valid = false;
                break;
            }
        }
        if (valid)
        {
            pathes.push_back(path);
        }
    }
    cout << "End filter short connections:\t" << pathes.size() << endl;
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes_buf4 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();

    for (auto path : pathes)
    {
        set<uint32_t> to_insert;
        for (int i = 0; i < path.size() - 1; i++)
        {
            to_insert.insert((*bubble_chain_begin_end_nodes_buf3)[path[i]][path[i + 1]].begin(), (*bubble_chain_begin_end_nodes_buf3)[path[i]][path[i + 1]].end());
        }
        if (to_insert.size() >= 2)
        {
            uint32_t beg = path[0];
            uint32_t end = path[path.size() - 1];
            if ((*bubble_chain_begin_end_nodes_buf4).find(beg) == (*bubble_chain_begin_end_nodes_buf4).end())
            {
                (*bubble_chain_begin_end_nodes_buf4)[beg] = map<uint32_t, set<uint32_t>>();
            }
            if ((*bubble_chain_begin_end_nodes_buf4).find(end ^ 1) == (*bubble_chain_begin_end_nodes_buf4).end())
            {
                (*bubble_chain_begin_end_nodes_buf4)[end ^ 1] = map<uint32_t, set<uint32_t>>();
            }
            (*bubble_chain_begin_end_nodes_buf4)[beg][end] = to_insert;
            (*bubble_chain_begin_end_nodes_buf4)[end ^ 1][beg ^ 1] = to_insert;
        }
    }
    delete bubble_chain_begin_end_nodes_buf3;
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes = new map<uint32_t, map<uint32_t, set<uint32_t>>>();

    for (auto beg : (*bubble_chain_begin_end_nodes_buf4))
    {
        (*bubble_chain_begin_end_nodes)[beg.first] = map<uint32_t, set<uint32_t>>();
        for (auto end : beg.second)
        {
            set<uint32_t> to_insert;
            for (auto n : end.second)
            {
                to_insert.insert(n >> 1);
            }
            uint64_t all_len = 0;
            // if(to_insert.size()>8 || (pure_incoming_num[beg.first]==0 && pure_outgoing_num[end.first] == 0)){
            for (auto n : to_insert)
            {
                all_len += g->seq[n].len;
            }
            // }
            // if(to_insert.size()>8 || all_len > 10000000 || (pure_incoming_num[beg.first]==0 && pure_outgoing_num[end.first] == 0) ){
            // if(to_insert.size()>8 || all_len > 10000000 || (pure_incoming_num[beg.first]==0 && pure_outgoing_num[end.first] == 0) ){
            // if(all_len > 10000000 ){
            if (all_len > 1000000)
            {
                (*bubble_chain_begin_end_nodes)[beg.first][end.first] = to_insert;
            }
        }
    }
    delete bubble_chain_begin_end_nodes_buf4;
    cout << "Filtered Nodes: " << endl;
    for (auto n : unaccessable)
    {
        cout << g->seq[n].name << ", ";
    }
    cout << endl;

    cout << "Not Filtered Nodes: " << endl;
    for (auto n : not_filtered)
    {
        cout << g->seq[n].name << ", ";
    }
    cout << endl;

    // for(auto a: *bubble_chain_begin_end_nodes){
    //     for(auto b: a.second){
    //         // all_nodes.insert(b.second.begin(),b.second.end());
    //         cout << g->seq[a.first>>1].name << ", " << g->seq[b.first>>1].name << endl;
    //     }
    // }

    uint32_t branch_counter = 0;
    set<uint32_t> all_nodes;
    for (auto a : *bubble_chain_begin_end_nodes)
    {
        for (auto b : a.second)
        {
            branch_counter++;
            all_nodes.insert(b.second.begin(), b.second.end());
        }
    }
    cout << "Nodes in original Graph: " << g->n_seq << endl;
    cout << "Total Nodes in Chain Graph: " << all_nodes.size() << endl;
    cout << "Branches: " << branch_counter << endl;

    for (auto bubble : bubbles)
    {
        delete bubble;
    }

    for (int i = 0; i < g->n_seq; i++)
    {
        free(connections_count[i]);
    }
    free(connections_count);
    return bubble_chain_begin_end_nodes;
    // cout << bubble_chain_end_begin.size() << endl;
    // cout << names.size() << endl;
    // cout << "finish get bubble chain" << endl;
}
typedef struct
{
    uint32_t begin, end;
    set<uint32_t> nodes;
    bool is_bubble;
    bool is_hap_only;
} bubble_chain;
struct NamedBubbleContig
{
    std::vector<std::string> hap1;
    std::vector<std::string> hap2;
};
typedef struct
{
    uint32_t begin, end;
    string hap_sequence;
    vector<uint32_t> hap_path;
    map<uint32_t, uint32_t> node_positions;
    vector<string> catched_contigs;
    bool is_valid;
} hap_chain_result_t;
vector<vector<bubble_chain>> bubble_chains;
std::vector<std::set<uint32_t>> components1;
vector<bubble_t *> bubbles_global;
int *node_type_global = NULL;
std::vector<std::set<uint32_t>> hap_chains_global;
vector<vector<hap_chain_result_t>> hap_results_debug;
vector<vector<uint32_t>> unvisited_nodes_global;
//..........................................................................

// void analyze_coverage(asg_t *g, int *node_type, int n_vtx)
// {
//     std::vector<int> coverage(n_vtx, 0);
//     for (int i = 0; i < n_vtx; i++)
//     {
//         coverage[i] = g->seq[i >> 1].coverage; //     }

//     auto [haploid_peak, diploid_peak] = estimate_diploid_coverage_levels(coverage.data(), n_vtx);

//     std::cout << "Estimated haploid coverage peak: " << haploid_peak << std::endl;
//     std::cout << "Estimated diploid coverage peak: " << diploid_peak << std::endl;
// }

void print_graph_edges(asg_t *graph)
{
    for (uint32_t v = 0; v < graph->n_seq * 2; ++v)
    {
        if (graph->seq[v >> 1].del)
            continue; 

        cout << "Node ID: " << v
             << ", Name: " << graph->seq[v >> 1].name
             << (v & 1 ? " (reverse)" : " (forward)") << endl;

        uint32_t out_degree = asg_arc_n(graph, v);
        asg_arc_t *out_arcs = asg_arc_a(graph, v);
        cout << "  Out-degree: " << out_degree << endl;
        for (uint32_t i = 0; i < out_degree; ++i)
        {
            uint32_t to = out_arcs[i].v;
            cout << "    → " << to
                 << " (" << graph->seq[to >> 1].name
                 << (to & 1 ? " [rev]" : " [fwd]") << ")" << endl;
        }

        uint32_t in_degree = asg_arc_n(graph, v ^ 1); // 入边 = 反向边的出边
        asg_arc_t *in_arcs = asg_arc_a(graph, v ^ 1);
        cout << "  In-degree: " << in_degree << endl;
        for (uint32_t i = 0; i < in_degree; ++i)
        {
            uint32_t from = in_arcs[i].v ^ 1;
            cout << "    ← " << from
                 << " (" << graph->seq[from >> 1].name
                 << (from & 1 ? " [rev]" : " [fwd]") << ")" << endl;
        }

        cout << endl;
    }
}

void deduplicate_outgoing_arcs_safe(asg_t *g)
{
    for (uint32_t v = 0; v < g->n_seq * 2; ++v)
    {
        uint32_t n = asg_arc_n(g, v);
        asg_arc_t *a = asg_arc_a(g, v);

        std::unordered_set<uint32_t> seen;

        for (uint32_t i = 0; i < n; ++i)
        {
            if (a[i].del)
                continue;

            uint32_t w = a[i].v;
            if (seen.find(w) != seen.end())
            {
                a[i].del = 1;             }
            else
            {
                seen.insert(w);
            }
        }
    }

    gfa_cleanup(g);
}

void asg_cut_arc(asg_t *g, uint32_t v, uint32_t w)
{
    asg_arc_t *av = asg_arc_a(g, v);
    int nv = asg_arc_n(g, v);
    for (int i = 0; i < nv; i++)
    {
        if (av[i].v == w)
        {
            av[i].del = 1;
        }
    }
}
vector<pair<uint32_t, uint32_t>> arcs_to_cut;
void build_bubbles(asg_t *g,
                   vector<bubble_t *> &bubbles,
                   vector<bubble_t *> bubble_by_ending_begining[],
                   int *node_type)
{
    bubbles.clear();
    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        bubble_by_ending_begining[i].clear();
        node_type[i] = 0;
    }

    uint32_t n_vtx = g->n_seq * 2;
    int delete_nums = 0;
    int bubbles_nums = 0;
    for (uint32_t i = 0; i < n_vtx; i++)
    {
        if (asg_arc_n(g, i) == 1 && asg_arc_n(g, i ^ 1) == 1)
        {
            asg_arc_t *out_arc = asg_arc_a(g, i);
            if (out_arc->del)
                continue;

            uint32_t v = out_arc->v;

            asg_arc_t *in_arc = asg_arc_a(g, i ^ 1);
            if (in_arc->del)
                continue;

            uint32_t v2 = in_arc->v;

            if ((v ^ 1) == v2)
            {
                asg_cut_arc(g, i, v);
                asg_cut_arc(g, i ^ 1, v2);
                asg_cut_arc(g, v, i);
                asg_cut_arc(g, v2, i ^ 1);
                delete_nums++;
                // cout << "delete: " << g->seq[i >> 1].name
                //      << " out: " << g->seq[v >> 1].name
                //      << " in: " << g->seq[v2 >> 1].name << endl;
                // gfa_cleanup(g);
            }
        }

        // if (i == 68826)
        // {
        //     cout << "test one :" << endl;
        //     cout << "Inspecting node i = " << i << " (" << g->seq[i >> 1].name
        //          << " [" << ((i & 1) ? "rev" : "fwd") << "])" << endl;

        //     asg_arc_t *out_arcs = asg_arc_a(g, i);
        //     uint32_t out_n = asg_arc_n(g, i);
        //     cout << "  Out-degree: " << out_n << endl;
        //     for (uint32_t k = 0; k < out_n; ++k)
        //     {
        //         uint32_t v = out_arcs[k].v;
        //         cout << "    → " << v << " (" << g->seq[v >> 1].name
        //              << " [" << ((v & 1) ? "rev" : "fwd") << "], del=" << out_arcs[k].del << ")" << endl;
        //     }

        //     asg_arc_t *in_arcs = asg_arc_a(g, i ^ 1);
        //     uint32_t in_n = asg_arc_n(g, i ^ 1);
        //     cout << "  In-degree: " << in_n << endl;
        //     for (uint32_t k = 0; k < in_n; ++k)
        //     {
        //         uint32_t u = in_arcs[k].v;
        //         cout << "    ← " << u << " (" << g->seq[u >> 1].name
        //              << " [" << ((u & 1) ? "rev" : "fwd") << "], del=" << in_arcs[k].del << ")" << endl;
        //     }
        // }

        bubble_t *result = detect_bubble(g, i);
        if (result != nullptr)
        {
            result->id = i;
            bubble_by_ending_begining[result->begNode].push_back(result);
            bubble_by_ending_begining[result->endNode].push_back(result);
            bubbles.push_back(result);
            bubbles_nums++;
            // std::cout << "Bubble detected with id = " << result->id
            //           << ", start node = " << g->seq[result->begNode / 2].name
            //           << ", end node = " << g->seq[result->endNode / 2].name << std::endl;
        }

        // if (i == 68826)
        // {
        //     cout << "test two :" << endl;
        //     cout << "Inspecting node i = " << i << " (" << g->seq[i >> 1].name
        //          << " [" << ((i & 1) ? "rev" : "fwd") << "])" << endl;

        //     asg_arc_t *out_arcs = asg_arc_a(g, i);
        //     uint32_t out_n = asg_arc_n(g, i);
        //     cout << "  Out-degree: " << out_n << endl;
        //     for (uint32_t k = 0; k < out_n; ++k)
        //     {
        //         uint32_t v = out_arcs[k].v;
        //         cout << "    → " << v << " (" << g->seq[v >> 1].name
        //              << " [" << ((v & 1) ? "rev" : "fwd") << "], del=" << out_arcs[k].del << ")" << endl;
        //     }

        //     asg_arc_t *in_arcs = asg_arc_a(g, i ^ 1);
        //     uint32_t in_n = asg_arc_n(g, i ^ 1);
        //     cout << "  In-degree: " << in_n << endl;
        //     for (uint32_t k = 0; k < in_n; ++k)
        //     {
        //         uint32_t u = in_arcs[k].v;
        //         cout << "    ← " << u << " (" << g->seq[u >> 1].name
        //              << " [" << ((u & 1) ? "rev" : "fwd") << "], del=" << in_arcs[k].del << ")" << endl;
        //     }
        // }
    }
    cout << "detect_bubble over. " << endl;
    for (auto &[src, tgt] : arcs_to_cut)
    {
        // cout << "[Before cut] Outgoing from " << g->seq[src >> 1].name << ": ";
        asg_arc_t *out_src = asg_arc_a(g, src);
        for (uint32_t i = 0; i < asg_arc_n(g, src); ++i)
        {
            if (!out_src[i].del)
            {
                cout << g->seq[out_src[i].v >> 1].name << " ";
            }
        }
        cout << endl;

        // cout << "[Before cut] Outgoing from " << g->seq[tgt >> 1].name << ": ";
        asg_arc_t *out_tgt = asg_arc_a(g, tgt ^ 1);
        for (uint32_t i = 0; i < asg_arc_n(g, tgt ^ 1); ++i)
        {
            if (!out_tgt[i].del)
            {
                // cout << g->seq[out_tgt[i].v >> 1].name << " ";
            }
        }
        // cout << endl;

        asg_cut_arc(g, src, tgt);

        // cout << "[After cut] Outgoing from " << g->seq[src >> 1].name << ": ";
        asg_arc_t *out_src1 = asg_arc_a(g, src);
        for (uint32_t i = 0; i < asg_arc_n(g, src); ++i)
        {
            if (!out_src1[i].del)
            {
                // cout << g->seq[out_src1[i].v >> 1].name << " ";
            }
        }
        // cout << endl;

        // cout << "[After cut] Outgoing from " << g->seq[tgt >> 1].name << ": ";
        asg_arc_t *out_tgt2 = asg_arc_a(g, tgt ^ 1);
        for (uint32_t i = 0; i < asg_arc_n(g, tgt ^ 1); ++i)
        {
            if (!out_tgt2[i].del)
            {
                // cout << g->seq[out_tgt2[i].v >> 1].name << " ";
            }
        }
        // cout << endl;
    }
    cout << "delete cycle : " << delete_nums << endl;
    cout << "bubble nums : " << bubbles_nums << endl;
    gfa_cleanup(g);

    // print_graph_edges(g);
    for (bubble_t *bubble : bubbles)
    {
        vector<asg_arc_t *> arc_stack;
        vector<uint32_t> node_stack, node_vi_stack;
        node_stack.push_back(bubble->begNode);
        node_vi_stack.push_back(0);
        // cout << "Now bubble: " << " begNode: " << g->seq[bubble->begNode / 2].name << " endNode: " << g->seq[bubble->endNode / 2].name << endl;
        int num = 0;
        while (!node_stack.empty())
        {
            num++;
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();

            if (u == bubble->endNode)
            {
                for (int ui = 0; ui < node_stack.size(); ui++)
                {
                    uint32_t nid = node_stack[ui];
                    if (ui == 0 || ui == node_stack.size() - 1)
                    {
                        if (node_type[nid] != 2)
                            node_type[nid] = 1;
                    }
                    else
                    {
                        node_type[nid] = 2;
                    }
                }
                // if (node_stack.size() == 3)
                // {
                //     uint32_t only_inner = node_stack[1];
                //     node_type[only_inner] = 4;
                // }

                bubble->paths.push_back(arc_stack);
                bubble->starting_arcs.insert(arc_stack.front()->v);
                bubble->ending_arcs.insert(node_stack[node_stack.size() - 2]);

                std::set<uint32_t> path_node_set;
                for (uint32_t nid : node_stack)
                {
                    path_node_set.insert(nid);
                }
                bubble->paths_nodes.push_back(path_node_set);

                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }

            // if (num > 10000)
            // {
            //     if (num % 1000 == 0)
            //         cout << "Now u is: " << g->seq[u / 2].name << std::endl;
            // }
            if (num > 2000000)
            {

                for (int ui = 0; ui < node_stack.size(); ui++)
                {
                    uint32_t nid = node_stack[ui];
                    node_type[nid] = 2;
                }
                cout << "Bubble path too long, break" << std::endl;
                node_type[bubble->begNode] = 1;
                node_type[bubble->endNode] = 1;
                break;
            }

            uint32_t num_out = asg_arc_n(g, u);
            if (vi < num_out)
            {
                asg_arc_t *out_arcs = asg_arc_a(g, u);
                node_vi_stack.back()++;
                arc_stack.push_back(out_arcs + vi);
                node_stack.push_back(out_arcs[vi].v);
                node_vi_stack.push_back(0);
            }
            else
            {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
        }
    }
}
unordered_set<uint32_t> complex_nodes;
void Get_bubble(asg_t *g,
                vector<bubble_t *> &bubbles,
                vector<bubble_t *> bubble_by_ending_begining[],
                int *node_type)
{

    build_bubbles(g, bubbles, bubble_by_ending_begining, node_type);

    // debug_filter_invalid_bubble_paths(g, bubbles);
    // cout << "node_type[i] = 3: " << endl;
    int type3_count = 0;
    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        if (node_type[i] != 1 && node_type[i] != 2)
        {
            uint32_t out_deg = asg_arc_n(g, i);    
            uint32_t in_deg = asg_arc_n(g, i ^ 1); 

            if (out_deg == 0 && in_deg != 0)
            {
                node_type[i] = 3;
                node_type[i ^ 1] = 3;
                // cout << g->seq[i >> 1].name << " ";
                type3_count++;
                // cout << node_type[i] << " ";
            }
            if (out_deg > 2 && in_deg > 2)
            {
                complex_nodes.insert(i);
            }
        }
        if (node_type[i] != 2)
        {
            uint32_t out_deg = asg_arc_n(g, i);    
            uint32_t in_deg = asg_arc_n(g, i ^ 1); 
            if (out_deg == 1 && in_deg == 1)
            {
                node_type[i] == 4;
            }
        }
    }
    cout << "type3_count: " << type3_count << endl;
}

vector<set<uint32_t>> find_connected_components_debug_plants(asg_t *g,
                                                             std::vector<std::set<uint32_t>> &hap_chains,
                                                             std::vector<std::set<uint32_t>> &filtered_components, int *node_type)
{
    vector<bool> visited_real_node(g->n_seq, false);
    vector<set<uint32_t>> components;

    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        uint32_t real_id = i >> 1;
        // if (visited_real_node[real_id]) continue;

        set<uint32_t> component;
        stack<uint32_t> s;
        s.push(i);

        while (!s.empty())
        {
            uint32_t u = s.top();
            s.pop();

            uint32_t uid = u >> 1;
            if (visited_real_node[uid])
                continue;
            visited_real_node[uid] = true;

            component.insert(uid * 2);    
            component.insert(uid * 2 + 1); 

            uint32_t out_deg = asg_arc_n(g, u);
            asg_arc_t *out_arcs = asg_arc_a(g, u);
            for (uint32_t j = 0; j < out_deg; j++)
            {
                uint32_t v = out_arcs[j].v;
                if (!visited_real_node[v >> 1])
                {
                    s.push(v);
                }
            }

  
            uint32_t in_deg = asg_arc_n(g, u ^ 1);
            asg_arc_t *in_arcs = asg_arc_a(g, u ^ 1);
            for (uint32_t j = 0; j < in_deg; j++)
            {
                uint32_t v = in_arcs[j].v;
                if (!visited_real_node[v >> 1])
                {
                    s.push(v);
                }
            }
        }


        uint32_t total_length = 0;
        for (auto node : component)
        {
            total_length += g->seq[node >> 1].len; 
        }
        hap_chains.push_back(component);
    }

    cout << endl;
    for (size_t i = 0; i < hap_chains.size(); ++i)
    {
        std::cout << "Haplotypic chain " << i << " contains nodes: ";
        for (auto node : hap_chains[i])
        {
            std::cout << g->seq[node >> 1].name << " ";
        }
        std::cout << std::endl;
    }

    for (size_t i = 0; i < components.size(); i++)
    {
        cout << "Component " << i << " contains nodes: ";
        vector<string> complex_node_names;
        int complex_count = 0;
        for (auto node : components[i])
        {
            cout << g->seq[node >> 1].name << " ";
            if (complex_nodes.count(node))
            {
                ++complex_count;
                complex_node_names.push_back(g->seq[node >> 1].name);
            }
        }
        cout << endl;
        if (complex_count > 0)
        {
            cout << "  >> Component " << i << " has " << complex_count << " complex nodes: ";
            if (complex_count > 10)
                complex_components.insert(i);
            for (const auto &name : complex_node_names)
                cout << name << " ";
            cout << endl;
        }
    }

    return components;
}

vector<set<uint32_t>> find_connected_components_debug(asg_t *g,
                                                      std::vector<std::set<uint32_t>> &hap_chains,
                                                      std::vector<std::set<uint32_t>> &filtered_components, int *node_type)
{
    vector<bool> visited_real_node(g->n_seq, false);
    vector<set<uint32_t>> components;

    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        uint32_t real_id = i >> 1;
        // if (visited_real_node[real_id]) continue;

        set<uint32_t> component;
        stack<uint32_t> s;
        s.push(i);

        while (!s.empty())
        {
            uint32_t u = s.top();
            s.pop();

            uint32_t uid = u >> 1;
            if (visited_real_node[uid])
                continue;
            visited_real_node[uid] = true;

            
            component.insert(uid * 2);     
            component.insert(uid * 2 + 1); 

 
            uint32_t out_deg = asg_arc_n(g, u);
            asg_arc_t *out_arcs = asg_arc_a(g, u);
            for (uint32_t j = 0; j < out_deg; j++)
            {
                uint32_t v = out_arcs[j].v;
                if (!visited_real_node[v >> 1])
                {
                    s.push(v);
                }
            }

            uint32_t in_deg = asg_arc_n(g, u ^ 1);
            asg_arc_t *in_arcs = asg_arc_a(g, u ^ 1);
            for (uint32_t j = 0; j < in_deg; j++)
            {
                uint32_t v = in_arcs[j].v;
                if (!visited_real_node[v >> 1])
                {
                    s.push(v);
                }
            }
        }

        uint32_t total_length = 0;
        for (auto node : component)
        {
            total_length += g->seq[node >> 1].len; 
        }

        if (component.size() >= 10)
        {
            components.push_back(component); 
        }
        else if (component.size() >= 4)
        {
            hap_chains.push_back(component); 
        }
    }

    for (auto &comp : components)
    {
        int total = comp.size();
        int type1_2_count = 0;
        int type4_count = 0;
        int begin_end = 0;
        double avg_degree = 0;
        int avg_degree_out = 0;
        int avg_degree_in = 0;
        for (auto node : comp)
        {
            uint32_t out_deg = asg_arc_n(g, node);
            avg_degree_out += out_deg;
            if (out_deg == 0)
            {
                begin_end++;
            }
            int type = node_type[node];
            if (type == 1 || type == 2)
            {
                type1_2_count++;
            }
            if (type == 4)
            {
                type4_count++;
            }
        }
        avg_degree = (double)avg_degree_out / total;
        double ratio = (double)type1_2_count / total; 
        double ratio2 = (double)type4_count / total;
        cout << "node: " << g->seq[*comp.begin() >> 1].name << " avg_degree: " << avg_degree << " type1_2_count : " << type1_2_count << " ratio " << ratio << " type4_count: " << type4_count << " ratio2: " << ratio2 << endl;
        // cout << "node: " << g->seq[*comp.begin() >> 1].name << " begin_end : " << begin_end << endl;
        if (ratio2 > 0.13)
        {
            hap_chains.push_back(comp);
        }

        else if (type1_2_count < 20 || ratio < 0.4 || (ratio < 0.61 && ratio2 >= 0.1))
        {
            if (avg_degree < 1.4)
                hap_chains.push_back(comp);
            else if (ratio2 >= 0.1 && ratio < 0.8)
            {
                hap_chains.push_back(comp);
            }
            else
            {
                filtered_components.push_back(comp);
            }
        }
        else if (total < 50 && ratio2 > 0.1 && ratio < 0.4)
        {
            hap_chains.push_back(comp);
        }
        // if (type1_2_count < 20 || ratio < 0.4  || type4_count > type1_2_count)
        // {
        //     cout<<"type1_2_count : "<<type1_2_count << " ratio "<<ratio<<" type4_count: "<<type4_count<<endl;
        //     hap_chains.push_back(comp);
        // }
        else
        {
            filtered_components.push_back(comp);
        }
    }

    components = std::move(filtered_components);
    cout << endl;

    for (size_t i = 0; i < hap_chains.size(); ++i)
    {
        std::cout << "Haplotypic chain " << i << " contains nodes: ";
        int count = 0;
        for (auto node : hap_chains[i])
        {
            std::cout << g->seq[node >> 1].name << " ";
            count++;
            if (count > 10)
            {
                break;
            }
        }
        std::cout << std::endl;
    }

    for (size_t i = 0; i < components.size(); i++)
    {
        cout << "Component " << i << " contains nodes: ";
        vector<string> complex_node_names;
        int complex_count = 0;
        int count = 0;
        for (auto node : components[i])
        {
            cout << g->seq[node >> 1].name << " ";
            if (complex_nodes.count(node))
            {
                ++complex_count;
                complex_node_names.push_back(g->seq[node >> 1].name);
                count++;
                if (count > 10)
                {
                    break;
                }
            }
        }
        cout << endl;
        if (complex_count > 0)
        {
            cout << "  >> Component " << i << " has " << complex_count << " complex nodes: ";
            if (complex_count > 10)
                complex_components.insert(i);
            for (const auto &name : complex_node_names)
                cout << name << " ";
            cout << endl;
        }
    }

    return components;
}

struct Short_branch_nodes
{
    uint32_t branch_begin;
    uint32_t branch_end;
    set<uint32_t> inside_nodes;
};
struct Type3Classification
{
    std::vector<Short_branch_nodes> short_branch_nodes;
    std::vector<uint32_t> chain_middle_nodes;                             // BRANCH
    std::vector<uint32_t> chain_end_nodes;                                // END
    std::set<std::pair<uint32_t, uint32_t>> potential_middle_connections; // TURN
};

void debug_classify_type3_nodes2(asg_t *g, int *node_type,
                                 Type3Classification &result,
                                 const std::unordered_set<uint32_t> &allowed_nodes)
{
    vector<uint32_t> chain_middle_nodes;
    vector<uint32_t> chain_end_nodes;

    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        if (allowed_nodes.count(i) == 0)
            continue;
        if (node_type[i] != 3)
            continue;

        // cout << "\n>>> Checking node ID: " << i << ", Name: " << g->seq[i >> 1].name << endl;
        uint32_t new_i = i;
        int num_incoming_arcs = asg_arc_n(g, i ^ 1);
        asg_arc_t *all_incoming_arcs = asg_arc_a(g, i ^ 1);
        if (num_incoming_arcs == 0)
        {
            // cout << "  Node " << i << " has no incoming arcs, skip." << endl;
            continue;
        }
        bool is_short_branch = true;
        bool is_branch = false;
        bool is_turn = false;
        bool is_END = true; 
        bool found_connect = false;
        uint32_t end_turn_node = (uint32_t)(-1);
        if (num_incoming_arcs == 2)
        {
            bool tip_found = false;
            for (uint32_t j = 0; j < num_incoming_arcs; j++)
            {
                uint32_t from = all_incoming_arcs[j].v;
                if (asg_arc_n(g, from) == 0)
                {
                    // cout << "node " << i << " has two incoming arcs, but one of them is a tip" << endl;
                    tip_found == true;
                    break;
                    ;
                }
            }
            if (tip_found)
            {
                continue;
            }
        }
        for (uint32_t j = 0; j < num_incoming_arcs; j++)
        {
            uint32_t from = all_incoming_arcs[j].v;
            // cout << "  Incoming arc from node: " << from << " (" << g->seq[from >> 1].name << ")" << endl;

            std::function<bool(uint32_t, std::set<uint32_t> &, int)> dfs;
            dfs = [&](uint32_t current, std::set<uint32_t> &visited, int depth) -> bool
            {
                if (depth >= 10)
                    return true;

                visited.insert(current);

                int next_deg = asg_arc_n(g, current);
                asg_arc_t *next_arcs = asg_arc_a(g, current);

                bool found_extend = false;
                for (int c = 0; c < next_deg; c++)
                {
                    uint32_t next = next_arcs[c].v;
                    if (visited.count(next))
                        continue;

                    // cout << "    DFS Step " << depth + 1 << ": From " << current << " to " << next
                    //      << " (" << g->seq[next >> 1].name << ")" << endl;

                    int fan_xiang = asg_arc_n(g, next ^ 1);
                    asg_arc_t *fan_xiang2 = asg_arc_a(g, next ^ 1);
                    for (int d = 0; d < fan_xiang; d++)
                    {
                        uint32_t found_fan = fan_xiang2[d].v;
                        if (asg_arc_n(g, found_fan) == 0)
                        {
                            // cout << " found only one node " << g->seq[found_fan >> 1].name << endl;
                            if (asg_arc_n(g, found_fan ^ 1) == 3)
                            {
                                // cout << "  Found turn node " << g->seq[found_fan >> 1].name << endl;
                                return found_extend;
                            }
                            found_connect = true;
                            result.potential_middle_connections.insert({i, found_fan ^ 1});
                        }
                    }

                    if (dfs(next, visited, depth + 1))
                    {
                        found_extend = true;
                        break; 
                    }
                }
                visited.erase(current);
                return found_extend;
            };

            // 从from开始做DFS
            std::set<uint32_t> visited;
            visited.insert(i ^ 1);
            visited.insert(i);
            bool extend5 = dfs(from ^ 1, visited, 0);

            if (extend5)
            {
                is_branch = true;
                is_END = false;
                // break; 
            }

            visited.clear();
            from = all_incoming_arcs[j].v;
            uint32_t current = from ^ 1;
            int step = 0;
            visited.insert(i ^ 1);
            visited.insert(i);
            if (num_incoming_arcs == 1)
            {
                current = asg_arc_a(g, i ^ 1)[0].v ^ 1;
                // cout << "current = " << current << " " << g->seq[current >> 1].name << endl;
            }
            while (step < 2)
            {
                int next_deg = asg_arc_n(g, current);
                if (next_deg == 0)
                {
                    // cout << "next_deg:" << asg_arc_n(g, current) << endl;
                    // cout << "  current node: " << current << " (" << g->seq[current >> 1].name << ")" << endl;
                    break;
                }

                asg_arc_t *next_arcs = asg_arc_a(g, current);
                bool advanced = false;

                for (int c = 0; c < next_deg; c++)
                {
                    uint32_t next = next_arcs[c].v;
                    if (visited.count(next))
                        continue;
                    current = next;
                    advanced = true;
                    break;
                }

                if (!advanced)
                {
                    // cout << "advanced 0" << endl;
                    break;
                }

                step++;

                if (step == 2)
                {
                    uint32_t turn_current = current ^ 1;
                    uint32_t turn_next_deg = asg_arc_n(g, turn_current);
                    asg_arc_t *turn_next_arcs = asg_arc_a(g, turn_current);
                    // cout << "  Turn node: " << turn_current << " (" << g->seq[turn_current >> 1].name << ")" << endl;
                    for (int c = 0; c < turn_next_deg; c++)
                    {
                        uint32_t turn_node = turn_next_arcs[c].v ^ 1;
                        if (node_type[turn_node] == 3 && turn_node != new_i)
                        {
                            uint32_t other_node = turn_next_arcs[1 - c].v ^ 1;
                            uint32_t a = g->seq[turn_node >> 1].len + g->seq[new_i >> 1].len;
                            uint32_t b = g->seq[other_node >> 1].len;
                            // cout << "      Checking TURN: " << turn_node << " (" << g->seq[turn_node >> 1].name << ")" << endl;
                            // cout << "      a: " << a << ", b: " << b
                            //      << ", diff: " << std::abs((int)a - (int)b)
                            //      << ", limit: " << 0.4 * std::max(a, b) << endl;
                            if (std::abs((int)a - (int)b) < 0.4 * std::max(a, b))
                            {
                                is_turn = true;
                                end_turn_node = turn_node;
                                if (asg_arc_n(g, turn_node ^ 1) == 2)
                                {
                                    is_turn = false;
                                    is_branch = true;
                                    // cout << "cant connecting" << endl;
                                }
                                if (asg_arc_n(g, turn_node) == 2)
                                {
                                    is_turn = false;
                                    is_branch = true;
                                    // cout << "cant connecting1" << endl;
                                }
                                // cout << "    ==> Found type3 node " << turn_node << " within 5 steps from " << new_i
                                //      << ", connecting: " << g->seq[new_i >> 1].name << " <--> " << g->seq[turn_node >> 1].name << endl;
                            }
                            else
                            {
                                is_branch = true;
                            }
                        }
                    }
                }
            }
        }

        if (is_branch && !is_turn)
        {
            // cout << "  ==> Node " << new_i << " is classified as BRANCH (middle)" << endl;
            result.chain_middle_nodes.push_back(i);
        }
        else if (is_turn)
        {
            // cout << "  ==> Node " << new_i << " is connecting " << end_turn_node << ", name: "
            //      << g->seq[i >> 1].name << " <--> " << g->seq[end_turn_node >> 1].name << endl;
            result.potential_middle_connections.insert({i, end_turn_node});
        }
        else if (is_END)
        {
            // cout << "  ==> Node " << new_i << " is classified as END" << endl;
            result.chain_end_nodes.push_back(i);
        }
    }
}
bool is_dead_end(asg_t *g, uint32_t start, int max_steps)
{
    std::set<uint32_t> visited;
    std::queue<std::pair<uint32_t, int>> q;
    q.push({start, 0});
    visited.insert(start);

    while (!q.empty())
    {
        auto [cur, steps] = q.front();
        q.pop();
        if (steps >= max_steps)
            return false; 

        int n = asg_arc_n(g, cur);
        asg_arc_t *a = asg_arc_a(g, cur);
        for (int i = 0; i < n; ++i)
        {
            uint32_t next = a[i].v;
            if (visited.count(next) == 0)
            {
                visited.insert(next);
                q.push({next, steps + 1});
            }
        }
    }

    return true; 
}

void add_edges_to_graph_new(asg_t *g, Type3Classification &type3Classification, int *node_type)
{

    std::set<std::pair<uint32_t, uint32_t>> to_remove1;
    for (const auto &p : type3Classification.potential_middle_connections)
    {
        std::pair<uint32_t, uint32_t> reversed = {p.second, p.first};
        if (p.first != p.second && type3Classification.potential_middle_connections.count(reversed))
        {
            if (p.first < p.second)
            {
                to_remove1.insert(reversed);
            }
        }
    }

    for (const auto &p : to_remove1)
    {
        type3Classification.potential_middle_connections.erase(p);
    }

    std::unordered_set<uint32_t> to_remove;

    std::unordered_map<uint32_t, uint32_t> pon_ends_map;

    for (const auto &b : type3Classification.chain_middle_nodes)
    {
        uint32_t v = b ^ 1;
        uint32_t nv = asg_arc_n(g, v);
        if (nv != 1)
        {

            continue; 
        }

        asg_arc_t *av = asg_arc_a(g, v);
        uint32_t w = av[0].v; 

        // if (node_type[w] != 1)
        // {
        //     continue;
        // }

        uint32_t w_rev = w ^ 1;
        uint32_t nw_rev = asg_arc_n(g, w_rev);
        if (nw_rev != 2)
        {

            continue;
        }

        asg_arc_t *aw_rev = asg_arc_a(g, w_rev);
        uint32_t v_rev = v ^ 1;
        uint32_t u = UINT32_MAX;
        for (uint32_t i = 0; i < nw_rev; ++i)
        {
            if (aw_rev[i].v != v_rev)
            {
                u = aw_rev[i].v;
                break;
            }
        }
        if (u == UINT32_MAX)
        {

            continue; 
        }

        uint32_t len_v = g->seq[v >> 1].len;
        uint32_t len_u = g->seq[u >> 1].len;

        double diff_ratio = std::abs((int64_t)len_v - (int64_t)len_u) / (double)std::max(len_v, len_u);
        // cout << "diff_ratio: " << diff_ratio << endl;
        int broken = 0;
        if (diff_ratio <= 0.15)
        {

            int step = 0;
            int type2_count = 0;
            std::set<uint32_t> visited;
            std::queue<uint32_t> q;
            q.push(u);
            visited.insert(u);
            to_remove.insert(b);
            while (!q.empty() && step < 100)
            {
                uint32_t cur = q.front();
                q.pop();
                ++step;

                // std::cout << "Step " << step << ": Visiting node "
                //           << g->seq[cur >> 1].name
                //           << " (node_type=" << node_type[cur] << ")" << std::endl;

                if (node_type[cur] == 2)
                    ++type2_count;

                uint32_t ncur = asg_arc_n(g, cur);
                asg_arc_t *acur = asg_arc_a(g, cur);
                if (asg_arc_n(g, cur ^ 1) == 2)
                {
                    for (int i = 0; i < asg_arc_n(g, cur ^ 1); i++)
                    {
                        uint32_t target = asg_arc_a(g, cur ^ 1)[i].v;
                        if (pon_ends_map.find(target) != pon_ends_map.end())
                        {
                            uint32_t matched_second = pon_ends_map[target];
                            broken = 1;
                            type3Classification.potential_middle_connections.insert({b, matched_second});
                            type3Classification.potential_middle_connections.insert({matched_second ^ 1, b ^ 1});
                        }
                    }
                }
                for (uint32_t i = 0; i < ncur; ++i)
                {
                    uint32_t next = acur[i].v;
                    if (visited.count(next) == 0)
                    {
                        visited.insert(next);
                        q.push(next);
                    }
                }
            }
            if (broken == 1)
            {
                break;
            }
            if (type2_count > 12)
            {
                asg_arc_t *au = asg_arc_a(g, u);
                for (uint32_t i = 0; i < asg_arc_n(g, u); ++i)
                {

                    type3Classification.potential_middle_connections.insert({b, au[i].v});
                    type3Classification.potential_middle_connections.insert({au[i].v ^ 1, b ^ 1});
                }
            }
            else
            {
                type3Classification.chain_end_nodes.push_back(v ^ 1);             }
        }
    }

    auto &nodes = type3Classification.chain_middle_nodes;
    nodes.erase(
        std::remove_if(nodes.begin(), nodes.end(),
                       [&](uint32_t n)
                       { return to_remove.count(n); }),
        nodes.end());

    for (const auto &c : type3Classification.chain_end_nodes)
    {
    }
    for (const auto &p : type3Classification.potential_middle_connections)
    {
        uint32_t u = p.first;
        uint32_t v = p.second;

        if ((u >> 1) >= g->n_seq || (v >> 1) >= g->n_seq)
        {
            std::cerr << "[Warning] Invalid node pair: " << u << ", " << v << std::endl;
            continue;
        }

        int32_t ov = 0, ow = 0;

        gfa_add_arc1(g, u, v, ov, ow, -1, 0);
        // gfa_add_arc1(g, v ^ 1, u ^ 1, ov, ow, -1, 0);

        // std::cout << "[Add Edge] " << g->seq[u >> 1].name << (u & 1 ? "-" : "+")
        //           << " -> " << g->seq[v >> 1].name << (v & 1 ? "-" : "+") << std::endl;
    }

    gfa_cleanup(g);
    for (const auto &v : type3Classification.chain_middle_nodes)
    {
        // cout << g->seq[(v ^ 1) >> 1].name << endl;
        int nv = asg_arc_n(g, v);
        asg_arc_t *av = asg_arc_a(g, v);

        for (int i = 0; i < nv; ++i)
        {
            uint32_t w = av[i].v;
            // cout << "w: " << g->seq[w >> 1].name << endl;
            int w_in_deg = 0;
            int nw_rev = asg_arc_n(g, w ^ 1);
            asg_arc_t *aw_rev = asg_arc_a(g, w ^ 1);
            for (int j = 0; j < nw_rev; ++j)
            {
                if ((aw_rev[j].v ^ 1) != v)
                    w_in_deg++;
            }

            if (w_in_deg > 0)
            {
                uint32_t current = w;
                set<uint32_t> visited;
                while (true)
                {
                    if (visited.count(current))
                        break;
                    visited.insert(current);
                    // cout << "current: " << g->seq[current >> 1].name << endl;
                    if (node_type[current] == 1 && node_type[current ^ 1] == 1)
                    {
                        uint32_t x = current;
                        // cout << "Found: " << g->seq[current >> 1].name << endl;
                        int nv_rev = asg_arc_n(g, v ^ 1);
                        asg_arc_t *av_rev = asg_arc_a(g, v ^ 1);
                        for (int k = 0; k < nv_rev; ++k)
                        {
                            uint32_t from = av_rev[k].v;
                            uint32_t len = g->seq[from >> 1].len;
                            uint32_t ol = av_rev[k].ol;

                            // cout << "[Redirect] " << g->seq[from >> 1].name << " -> "
                            //      << g->seq[x >> 1].name << endl;
                            // cout << "from: " << from << " x: " << x << endl;
                            gfa_add_arc1(g, from, x, 0, 0, -1, 0);
                            gfa_add_arc1(g, x ^ 1, from, 0, 0, -1, 0);
                            // add_edge_to_graph(g, from, x, len, ol);
                            // add_edge_to_graph(g, x ^ 1, from ^ 1, len, ol);
                        }
                        break;
                    }

                    int n_next = asg_arc_n(g, current);
                    asg_arc_t *a_next = asg_arc_a(g, current);
                    if (n_next == 0)
                        break;

                    current = a_next[0].v;
                }
            }
        }
    }

    for (const auto &v : type3Classification.chain_end_nodes)
    {
        int nv = asg_arc_n(g, v ^ 1);
        asg_arc_t *av = asg_arc_a(g, v);
        if (nv == 1)
        {
            uint32_t w = av[0].v;
            uint32_t w1 = w ^ 1;
            int nw = asg_arc_n(g, w1);
            asg_arc_t *aw = asg_arc_a(g, w1);
            uint32_t onther_end_node = 0;

            if (nw == 2)
            {
                bool real_end_node = false;
                for (int i = 0; i < nw; ++i)
                {
                    uint32_t reward = aw[i].v ^ 1;
                    if (reward == v)
                    {
                        real_end_node = true;
                        continue;
                    }
                    else
                    {
                        onther_end_node = reward;
                    }
                }

                if (real_end_node == true)
                {
                    // cout << "is real end: " << v << " " << g->seq[v >> 1].name << endl;
                    continue;
                }
            }
        }

        if (nv == 2)
        {
            for (int i = 0; i < nv; ++i)
            {
                uint32_t rv = av[i].v ^ 1;
                int nw = asg_arc_n(g, rv);
                asg_arc_t *aw = asg_arc_a(g, rv);

                if (nw == 2)
                {
                    for (int i = 0; i < nw; ++i)
                    {
                        uint32_t reward = aw[i].v ^ 1;
                        if (reward == v)
                            continue;

                        if (is_dead_end(g, reward, 3))
                        {
                            // cout << "is real end: " << v << " " << g->seq[v >> 1].name << endl;
                            break;
                        }
                    }
                }
            }
        }

    }

    gfa_cleanup(g);
}
bubble_t *detect_bubble3(asg_t *g, uint32_t source)
{
    uint32_t num_source_outgoing_arcs = asg_arc_n(g, source);
    uint32_t n_vtx = g->n_seq * 2;
    uint32_t num_seen = 1;
    if (num_source_outgoing_arcs < 2)
    {
        return nullptr;
    }

    bool seen[n_vtx], visited[n_vtx];
    for (int u = 0; u < n_vtx; u++)
    {
        seen[u] = false;
        visited[u] = false;
    }
    vector<uint32_t> frontier;

    frontier.push_back(source);
    seen[source] = true;

    while (!frontier.empty())
    {
        uint32_t v = frontier.back();
        frontier.pop_back();
        // cout << v << "\t" << g->seq[v / 2].name << endl; // topological order
        // cout<<"seen[v]: "<<seen[v]<<" visited[v]: "<<visited[v]<<endl;
        // assert(seen[v]);
        if (visited[v])
        {
            return nullptr;
        }
        // assert(!visited[v]);
        visited[v] = true;
        seen[v] = false;
        num_seen--;
        uint32_t num_outgoing_arcs = asg_arc_n(g, v);

        if (num_outgoing_arcs == 0)
        {
            return nullptr;
        }

        // add children to frontier if they have no other unvisited parents
        // uint32_t num_outgoing_arcs = asg_arc_n(g, v);  // 2
        asg_arc_t *outgoing_arcs = asg_arc_a(g, v); // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
        set<uint32_t> edge_set;
        for (int vi = 0; vi < num_outgoing_arcs; vi++)
        {
            uint32_t u = outgoing_arcs[vi].v;
            if (edge_set.insert(u).second)
            {
                if (u == v || u == (v ^ 1) || visited[u ^ 1] || u == source)
                {
                    return nullptr;
                }

                // assert(!visited[u]);
                if (!seen[u])
                {
                    num_seen++;
                }
                seen[u] = true;
                // assert(asg_arc_n(g, u ^ 1) >= 1);
                bool has_unvisited_parents = false;
                uint32_t num_incoming_arcs = asg_arc_n(g, u ^ 1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u ^ 1);
                for (int viuj = 0; viuj < num_incoming_arcs; viuj++)
                {
                    uint32_t p = (incoming_arcs_complement[viuj].v) ^ 1; // parent
                    if (!visited[p])
                        has_unvisited_parents = true;
                }

                if (!has_unvisited_parents)
                {
                    frontier.push_back(u);
                }
            }
        }

        if (frontier.size() == 1 && num_seen == 1 && seen[frontier.front()])
        {
            uint32_t t = frontier.back();
            frontier.pop_back();
            uint32_t num_outgoing_arcs = asg_arc_n(g, t);
            if (num_outgoing_arcs > 0)
            {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, t);
                for (int vi = 0; vi < num_outgoing_arcs; vi++)
                {
                    if (outgoing_arcs[vi].v == source)
                    {
                        return nullptr;
                    }
                }
            }
            // bool has_direct_path = false;
            // asg_arc_t *out_arcs = asg_arc_a(g, source);
            // for (int i = 0; i < asg_arc_n(g, source); ++i)
            // {
            //     if (out_arcs[i].v == t)
            //     {
            //         has_direct_path = true;
            //         cout << "has_direct_path: " << g->seq[source >> 1].name << " " << g->seq[t >> 1].name << endl;
            //         break;
            //     }
            // }

            // if (has_direct_path)
            // {
            //     uint32_t source_out_deg = asg_arc_n(g, source);
            //     uint32_t t_in_deg = asg_arc_n(g, t ^ 1);
            //     if (source_out_deg == 2 && t_in_deg == 2)
            //     {
            //         cout << "Cutting direct arc " << g->seq[source >> 1].name << " -> " << g->seq[t >> 1].name
            //              << endl;
            //         // asg_cut_arc(g, source, t);
            //         arcs_to_cut.emplace_back(source, t);
            //         return nullptr;
            //     }
            //     else
            //     {
            //         cout << "Not cutting arc: source_out_deg=" << source_out_deg << ", t_in_deg=" << t_in_deg << endl;
            //     }
            // }
            bubble_t *result = new bubble_t();
            result->begNode = source;
            result->endNode = t;
            return result;
        }
    }

    return nullptr;
}

void build_bubbles2(asg_t *g,
                    vector<bubble_t *> &bubbles,
                    vector<bubble_t *> bubble_by_ending_begining[],
                    int *node_type)
{
    // bubbles.clear();
    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        bubble_by_ending_begining[i].clear();
        node_type[i] = 0;
    }

    uint32_t n_vtx = g->n_seq * 2;
    for (uint32_t i = 0; i < n_vtx; i++)
    {
        // cout << i << " " << g->seq[i / 2].name << " " << endl;
        bubble_t *result = detect_bubble3(g, i);
        if (result != nullptr)
        {
            result->id = i;
            bubble_by_ending_begining[result->begNode].push_back(result);
            bubble_by_ending_begining[result->endNode].push_back(result);
            bubbles.push_back(result);
            // std::cout << "Bubble detected with id = " << result->id
            //           << ", start node = " << g->seq[result->begNode / 2].name
            //           << ", end node = " << g->seq[result->endNode / 2].name << std::endl;
        }
    }

    cout << "node_type[i] = 3: " << endl;
    for (uint32_t i = 0; i < g->n_seq * 2; i++)
    {
        if (node_type[i] != 2)
        {
            uint32_t out_deg = asg_arc_n(g, i);    
            uint32_t in_deg = asg_arc_n(g, i ^ 1); 

            if (out_deg == 0 && in_deg != 0)
            {
                node_type[i] = 3;
                node_type[i ^ 1] = 3;
                cout << g->seq[i >> 1].name << " ";

                // cout << node_type[i] << " ";
            }
            if (out_deg > 2 && in_deg > 2)
            {
                complex_nodes.insert(i);
            }
        }
        // if (node_type[i] != 2)
        // {
        //     uint32_t out_deg = asg_arc_n(g, i);    
        //     uint32_t in_deg = asg_arc_n(g, i ^ 1); 
        //     if (out_deg == 1 && in_deg == 1)
        //     {
        //         node_type[i] == 4;
        //     }
        // }
    }

    for (int b = 0; b < bubbles.size(); b++)
    {
        bubble_t *bubble = bubbles[b];
        uint32_t bubble_beginning = bubble->begNode;
        uint32_t bubble_end = bubble->endNode;

        // cout << "start get bubble paths from " << g->seq[bubble_beginning / 2].name << " to " << g->seq[bubble_end / 2].name << endl;

        vector<asg_arc_t *> arc_stack;
        vector<uint32_t> node_stack; // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_beginning);
        node_vi_stack.push_back(0);
        int stack_count = 0;
        uint32_t lastNode = bubble_beginning;
        int num = 0;
        while (!node_stack.empty())
        {
            num++;
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            if (u == bubble_end)
            {
                for (int ui = 0; ui < node_stack.size(); ui++)
                {
                    // cout << g->seq[node_stack[ui]/2].name << " ";
                    if (ui == 0 || ui == node_stack.size() - 1)
                    {
                        if (node_type[node_stack[ui]] != 2)
                        {
                            node_type[node_stack[ui]] = BUBBLE_END_BEGIN;
                        }
                    }
                    else
                    {
                        node_type[node_stack[ui]] = BUBBLE_INSIDE;
                    }
                }
                stack_count++;
                // cout << endl;

                bubble->paths.push_back(arc_stack);
                // cout << "Bubble " << b << " has " << bubble->paths.size() << " paths" << endl;
                // for (int p=0; p < bubble->paths.size(); p++) {
                //     cout << "  Path " << p << " length: " << bubble->paths[p].size() << endl;
                // }
                bubble->starting_arcs.insert(arc_stack.front()->v);
                bubble->ending_arcs.insert(lastNode);

                std::set<uint32_t> path_node_set;
                for (uint32_t nid : node_stack)
                {
                    path_node_set.insert(nid);
                }
                bubble->paths_nodes.push_back(path_node_set);
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }
            if (num > 2000000)
            {

                for (int ui = 0; ui < node_stack.size(); ui++)
                {
                    uint32_t nid = node_stack[ui];
                    node_type[nid] = 2;
                }
                cout << "Bubble path too long, break" << std::endl;
                node_type[bubble->begNode] = 1;
                node_type[bubble->endNode] = 1;
                break;
            }
            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            if (vi < num_outgoing_arcs)
            {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
                uint32_t v = outgoing_arcs[vi].v;
                node_vi_stack.back()++;
                arc_stack.push_back(outgoing_arcs + vi);
                node_stack.push_back(v);
                node_vi_stack.push_back(0);
            }
            else
            {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
            lastNode = u;
        }
        // cout << stack_count << endl;
        // cout << "finish get bubble paths from " << bubble_beginning << " to " << bubble_end << endl;
    }
}

void analyze_bubble_node_degrees_debug(
    asg_t *g,
    int *node_type,
    uint32_t n_vtx,
    vector<bubble_t *> bubble_by_ending_begining[],
    set<uint32_t> &bubble_chain_end_begin,
    map<uint32_t, uint32_t> &pure_outgoing_num,
    map<uint32_t, uint32_t> &pure_incoming_num,
    set<uint32_t> &out_only,
    const std::vector<std::set<uint32_t>> &components) {
    for (const auto &component : components)
    {
        for (uint32_t i : component) 
        {
            if (node_type[i] != 2)
            {
                int num_outgoing_arcs = asg_arc_n(g, i);
                int num_incoming_arcs = asg_arc_n(g, i ^ 1);

                asg_arc_t *all_outgoing_arcs = asg_arc_a(g, i);
                asg_arc_t *all_incoming_arcs = asg_arc_a(g, i ^ 1);

                set<uint32_t> valid_outgoing_arcs, valid_incoming_arcs;
                for (int c = 0; c < num_outgoing_arcs; c++)
                    valid_outgoing_arcs.insert(all_outgoing_arcs[c].v);
                for (int c = 0; c < num_incoming_arcs; c++)
                    valid_incoming_arcs.insert(all_incoming_arcs[c].v ^ 1);

                num_incoming_arcs = valid_incoming_arcs.size();
                num_outgoing_arcs = valid_outgoing_arcs.size();

                if (node_type[i] == 1)
                {
                    vector<bubble_t *> &bubbles = bubble_by_ending_begining[i];
                    set<uint32_t> outgoing_arcs_set, incoming_arcs_set;

                    sort(bubbles.begin(), bubbles.end(), [](const auto &lhs, const auto &rhs)
                         { return lhs->starting_arcs.size() > rhs->starting_arcs.size(); });

                    for (auto *bubble : bubbles)
                    {
                        if (bubble->begNode == i)
                        {
                            bool flag = false;
                            for (auto c : bubble->starting_arcs)
                            {
                                if (outgoing_arcs_set.insert(c).second)
                                {
                                    num_outgoing_arcs--;
                                    flag = true;
                                }
                            }
                            if (flag)
                                num_outgoing_arcs++;
                        }
                    }

                    sort(bubbles.begin(), bubbles.end(), [](const auto &lhs, const auto &rhs)
                         { return lhs->ending_arcs.size() > rhs->ending_arcs.size(); });

                    for (auto *bubble : bubbles)
                    {
                        if (bubble->endNode == i)
                        {
                            bool flag = false;
                            for (auto c : bubble->ending_arcs)
                            {
                                if (incoming_arcs_set.insert(c).second)
                                {
                                    num_incoming_arcs--;
                                    flag = true;
                                }
                            }
                            if (flag)
                                num_incoming_arcs++;
                        }
                    }
                }

                pure_outgoing_num[i] = num_outgoing_arcs;
                pure_incoming_num[i] = num_incoming_arcs;
                if (asg_arc_n(g, i ^ 1) == 1 && asg_arc_n(g, i) == 2 && num_outgoing_arcs == 1 && num_incoming_arcs == 1)
                {
                    // cout << "Debug: " << g->seq[i / 2].name << ": incoming " << num_incoming_arcs << "; outgoing " << num_outgoing_arcs
                    //      << "; bubbles: " << bubble_by_ending_begining[i].size() << endl;
                    bubble_chain_end_begin.insert(i);
                    bubble_chain_end_begin.insert(i ^ 1);
                }
                if (num_incoming_arcs == 0 && num_outgoing_arcs > 0)
                    out_only.insert(i);

                if (!(num_incoming_arcs == 1 && num_outgoing_arcs == 1))
                {
                    // cout << g->seq[i / 2].name << ": incoming " << num_incoming_arcs << "; outgoing " << num_outgoing_arcs
                    //      << "; bubbles: " << bubble_by_ending_begining[i].size() << endl;
                    bubble_chain_end_begin.insert(i);
                }
                if ((num_incoming_arcs == 1 && num_outgoing_arcs == 1) && bubble_by_ending_begining[i ^ 1].size() == 1 && bubble_by_ending_begining[i ^ 1].size() == 0)
                {
                    // cout << "fake bubble: " << g->seq[i / 2].name << ": incoming " << num_incoming_arcs << "; outgoing " << num_outgoing_arcs
                    //      << "; bubbles: " << bubble_by_ending_begining[i].size() << endl;
                    bubble_chain_end_begin.insert(i);
                }
            }
        }
    }
}

void write_gfa_file(const string &filename, asg_t *g, const vector<vector<bubble_chain>> &bubble_chains, bool only_bubble)
{
    ofstream out(filename);
    unordered_set<uint32_t> printed_nodes;

    for (const auto &component : bubble_chains)
    {
        for (const auto &bc : component)
        {
            if (only_bubble && !bc.is_bubble)
                continue;
            if (!only_bubble && !(bc.is_bubble || bc.is_hap_only))
                continue;

            for (uint32_t node : {bc.begin, bc.end})
            {
                if (printed_nodes.count(node))
                    continue;
                printed_nodes.insert(node);
                out << "S\t" << g->seq[node >> 1].name << "\t*\tLN:i:" << g->seq[node >> 1].len << "\n";
            }
        }
    }

    for (const auto &component : bubble_chains)
    {
        for (const auto &bc : component)
        {
            if (only_bubble && !bc.is_bubble)
                continue;
            if (!only_bubble && !(bc.is_bubble || bc.is_hap_only))
                continue;

            out << "L\t" << g->seq[bc.begin >> 1].name << "\t" << (bc.begin % 2 == 0 ? "+" : "-") << "\t"
                << g->seq[bc.end >> 1].name << "\t" << (bc.end % 2 == 0 ? "+" : "-") << "\t0M\n";
        }
    }

    out.close();
}
void extract_bubble_chains_by_dfs(asg_t *g,
                                  const set<uint32_t> &bubble_chain_begin_end,
                                  int *node_type,
                                  const vector<set<uint32_t>> &components,
                                  vector<vector<bubble_chain>> &bubble_chains, string output_directory)
{
    bubble_chains.resize(components.size());
    for (auto bubble_chain_begin : bubble_chain_begin_end)
    {
        map<uint32_t, set<uint32_t>> current_map;
        vector<asg_arc_t *> arc_stack;
        vector<uint32_t> node_stack; // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_chain_begin);
        node_vi_stack.push_back(0);
        bool visited[g->n_seq * 2];
        for (int u = 0; u < g->n_seq * 2; u++)
        {
            visited[u] = false;
        }

        while (!node_stack.empty())
        {
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
            visited[u] = true;
            int flag = 0;
            if (u != bubble_chain_begin && find(bubble_chain_begin_end.begin(), bubble_chain_begin_end.end(), u) != bubble_chain_begin_end.end())
            {
                flag = 1;
            }
            else
            {
                for (auto b : current_map)
                {
                    if (flag)
                    {
                        break;
                    }
                    for (int i = 0; i < num_outgoing_arcs; i++)
                    {
                        if (b.second.find(outgoing_arcs[i].v) != b.second.end())
                        {
                            u = b.first;
                            flag = 2;
                            break;
                        }
                    }
                }
            }

            if (flag == 1)
            {
                if (current_map.find(u) == current_map.end())
                {
                    set<uint32_t> new_set;
                    current_map.insert({u, new_set});
                }
                if (arc_stack.size() > 0)
                {
                    for (int a = 0; a < arc_stack.size() - 1; a++)
                    {
                        current_map[u].insert(arc_stack[a]->v);
                    }
                }
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }
            else if (flag == 2)
            {
                if (current_map.find(u) == current_map.end())
                {
                    set<uint32_t> new_set;
                    current_map.insert({u, new_set});
                }
                if (arc_stack.size() > 0)
                {
                    for (int a = 0; a < arc_stack.size(); a++)
                    {
                        current_map[u].insert(arc_stack[a]->v);
                    }
                }
            }

            if (vi < num_outgoing_arcs)
            {
                while (vi < num_outgoing_arcs)
                {
                    uint32_t v = outgoing_arcs[vi].v;
                    node_vi_stack.back()++;
                    if (!visited[v])
                    {
                        arc_stack.push_back(outgoing_arcs + vi);
                        node_stack.push_back(v);
                        node_vi_stack.push_back(0);
                        break;
                    }
                    else
                    {
                        vi++;
                    }
                }
            }
            if (vi >= num_outgoing_arcs)
            {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
        }
        {
            set<uint32_t> path_nodes;
            int bubble_count = 0;
            int hap_count = 0;
            for (auto map : current_map)
            {
                uint32_t end = map.first;
                for (auto node : current_map[end])
                {
                    path_nodes.insert(node);
                    if (node_type[node] == 1)
                        bubble_count++;
                    if (node_type[node] == 4)
                        hap_count++;
                }

                bubble_chain bc;
                bc.begin = bubble_chain_begin;
                bc.end = end;
                bc.nodes = path_nodes;
                bc.is_hap_only = bubble_count < 3;
                if (bubble_count > 3 && bubble_count < 10 && hap_count > 3)
                    bc.is_hap_only = true;
                // bc.is_hap_only = hap1_count > 1 && (hap1_count * 1.0 / total_node_count >= 0.9);
                bc.is_bubble = !bc.is_hap_only;
                // cout << " chain: " << g->seq[bubble_chain_begin >> 1].name << "  " << g->seq[end >> 1].name << " bubble_count: " << bubble_count << endl;
                for (size_t i = 0; i < components.size(); ++i)
                {
                    if (components[i].count(bubble_chain_begin))
                    {
                        bubble_chains[i].push_back(bc);
                        break;
                    }
                }
            }
        }

        // bubble_chain_begin_end_nodes->insert({bubble_chain_begin, current_map});
    }

    for (auto &component_chains : bubble_chains)
    {
        for (auto &bc : component_chains)
        {
            bc.nodes.insert(bc.begin);
            bc.nodes.insert(bc.end);
        }
    }

    system((string("rm -r ") + output_directory).c_str());
    mkdir(output_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    write_gfa_file(output_directory + "/main_graph_bubble_only.gfa", g, bubble_chains, true);

    write_gfa_file(output_directory + "/main_graph.gfa", g, bubble_chains, false);

    ofstream node_colors(output_directory + string("/") + string("main_graph_colors.csv"));
    node_colors << "NodeName,Color" << endl; // CSV 表头

    unordered_set<uint32_t> printed_nodes;
    unordered_set<uint32_t> bubble_endpoints;
    for (const auto &component : bubble_chains)
    {
        for (const auto &bc : component)
        {
            if (!bc.is_bubble)
                continue;
            printed_nodes.insert(bc.begin);
            printed_nodes.insert(bc.end);

            if (bc.is_bubble)
            {
                bubble_endpoints.insert(bc.begin);
                bubble_endpoints.insert(bc.end);
            }
        }
    }
    ofstream out(output_directory + string("/") + string("main_graph_colors.csv"));
    if (!out)
    {
        cerr << "ERROR: cannot open color file " << output_directory << " for writing\n";
        return;
    }

    out << "NodeName,Color\n";
    vector<uint32_t> nodes;
    nodes.reserve(printed_nodes.size());
    for (auto n : printed_nodes)
        nodes.push_back(n);
    sort(nodes.begin(), nodes.end());

    for (auto node : nodes)
    {
        string name = g->seq[node >> 1].name;
        string color = (bubble_endpoints.count(node) ? "red" : "black");
        out << name << "," << color << "\n";
    }
    out.close();

    cout << "finish get bubble chain" << endl;
}

void collect_bubble_chains_to_map(
    const std::vector<std::vector<bubble_chain>> &bubble_chains,
    std::map<uint32_t, std::map<uint32_t, std::set<uint32_t>>> *bubble_chain_begin_end_nodes_buf)
{
    for (const auto &component_chains : bubble_chains)
    {
        for (const auto &bc : component_chains)
        {
            if (!bc.is_bubble)
                continue;
            (*bubble_chain_begin_end_nodes_buf)[bc.begin][bc.end] = bc.nodes;
        }
    }
}
vector<uint32_t> bfs_find_path(asg_t *g, const set<uint32_t> &node_set,
                               uint32_t start, unordered_set<uint32_t> &global_visited);
void find_all_paths_in_hap_chains_BFS(asg_t *g, const vector<set<uint32_t>> &hap_chains, int *node_type)
{
    for (size_t chain_idx = 0; chain_idx < hap_chains.size(); ++chain_idx)
    {
        const auto &node_set = hap_chains[chain_idx];

        cout << "\n=== Chain " << chain_idx << " (" << node_set.size() << " nodes) ===" << endl;

        unordered_set<uint32_t> global_visited;
        vector<vector<uint32_t>> all_paths;     

        for (uint32_t start : node_set)
        {
            if (global_visited.count(start))
                continue;

            vector<uint32_t> current_path = bfs_find_path(g, node_set, start, global_visited);

            if (!current_path.empty() && current_path.size() > 1)
            {
                all_paths.push_back(current_path);

                for (uint32_t node : current_path)
                {
                    global_visited.insert(node);
                    global_visited.insert(node ^ 1); 
                }
            }
        }

        vector<uint32_t> longest_path;
        if (!all_paths.empty())
        {
            sort(all_paths.begin(), all_paths.end(),
                 [](const vector<uint32_t> &a, const vector<uint32_t> &b)
                 {
                     return a.size() > b.size(); 
                 });

            longest_path = all_paths[0];

            cout << "\nLongest path in chain " << chain_idx << " (length = "
                 << longest_path.size() << "):" << endl;
            for (size_t i = 0; i < longest_path.size(); ++i)
            {
                uint32_t node = longest_path[i];
                cout << "  [" << i + 1 << "] "
                     << g->seq[node >> 1].name
                     << (node & 1 ? "-" : "+");

                if (node_type[node] >= 0)
                    cout << " (type: " << node_type[node] << ")";
                cout << endl;
            }
        }

        if (!longest_path.empty())
        {
            hap_chain_result_t result;
            result.begin = longest_path.front();
            result.end = longest_path.back();
            result.hap_path = longest_path;
        }

        vector<uint32_t> unvisited_nodes;
        for (uint32_t node : node_set)
        {
            if (!global_visited.count(node))
            {
                unvisited_nodes.push_back(node);
                cout << "Unvisited node: " << g->seq[node >> 1].name
                     << (node & 1 ? "-" : "+")
                     << " (type: " << node_type[node] << ")" << endl;
            }
        }

        cout << "\nChain " << chain_idx << " summary:"
             << "\n  Total nodes: " << node_set.size()
             << "\n  Found paths: " << all_paths.size()
             << "\n  Longest path length: " << (longest_path.empty() ? 0 : longest_path.size())
             << "\n  Unvisited nodes: " << unvisited_nodes.size() << endl;

        unvisited_nodes_global.push_back(unvisited_nodes);
    }
}

vector<uint32_t> bfs_find_path(asg_t *g, const set<uint32_t> &node_set,
                               uint32_t start, unordered_set<uint32_t> &global_visited)
{
    unordered_map<uint32_t, uint32_t> parent; 
    unordered_set<uint32_t> visited;          
    queue<uint32_t> q;
    vector<uint32_t> end_nodes; 

    visited.insert(start);
    q.push(start);
    parent[start] = UINT32_MAX;

    while (!q.empty())
    {
        uint32_t u = q.front();
        q.pop();

        bool has_next = false;
        uint32_t deg = asg_arc_n(g, u);
        asg_arc_t *arcs = asg_arc_a(g, u);

        for (uint32_t i = 0; i < deg; ++i)
        {
            if (arcs[i].del)
                continue;

            uint32_t v = arcs[i].v;

            if (node_set.count(v) && !global_visited.count(v) && !visited.count(v))
            {
                visited.insert(v);
                parent[v] = u;
                q.push(v);
                has_next = true;
            }
        }

        if (!has_next)
        {
            end_nodes.push_back(u);
        }
    }

    if (end_nodes.empty())
        return {};

    uint32_t end_node = end_nodes[0];

    vector<uint32_t> path;
    uint32_t current = end_node;
    while (current != UINT32_MAX)
    {
        path.push_back(current);
        current = parent[current];
    }

    reverse(path.begin(), path.end());

    return path;
}

int check_connection(asg_t *g, uint32_t from_node, std::unordered_map<uint32_t, int> begin2chains_id, std::unordered_map<uint32_t, int> end2chains_id, unordered_set<uint32_t> &used)
{

    uint32_t deg = asg_arc_n(g, from_node);
    asg_arc_t *arcs = asg_arc_a(g, from_node);

    for (uint32_t i = 0; i < deg; ++i)
    {
        if (arcs[i].del)
            continue;

        if (begin2chains_id.find(arcs[i].v >> 1) != begin2chains_id.end())
        {
            int chain_id = begin2chains_id.at(arcs[i].v >> 1);
            if (used.count(chain_id))
            {
                cout << "chain_id is used:" << chain_id << endl;
            }
            else
            {
                return 1;
            }
        }

        if (end2chains_id.find(arcs[i].v >> 1) != end2chains_id.end())
        {
            int chain_id = end2chains_id.at(arcs[i].v >> 1);
            if (used.count(chain_id))
            {
                cout << "chain_id is used:" << chain_id << endl;
            }
            else
            {
                return 2;
            }
        }
    }

    uint32_t deg1 = asg_arc_n(g, from_node ^ 1);
    asg_arc_t *arcs1 = asg_arc_a(g, from_node ^ 1);

    for (uint32_t i = 0; i < deg1; ++i)
    {
        if (arcs1[i].del)
            continue;

        if (begin2chains_id.find(arcs1[i].v >> 1) != begin2chains_id.end())
        {
            int chain_id = begin2chains_id.at(arcs1[i].v >> 1);
            if (used.count(chain_id))
            {
                cout << "chain_id is used:" << chain_id << endl;
            }
            else
            {
                return 1;
            }
        }

        if (end2chains_id.find(arcs1[i].v >> 1) != end2chains_id.end())
        {
            int chain_id = end2chains_id.at(arcs1[i].v >> 1);
            if (used.count(chain_id))
            {
                cout << "chain_id is used:" << chain_id << endl;
            }
            else
            {
                return 2;
            }
        }
    }
    return 0;
}

string reverse_complement(string unitig)
{
    stringstream result;
    int c = 0;
    while (c < unitig.size())
    {
        char current = unitig[c];
        if (current == 'A')
        {
            result << 'T';
        }
        else if (current == 'T')
        {
            result << 'A';
        }
        else if (current == 'C')
        {
            result << 'G';
        }
        else if (current == 'G')
        {
            result << 'C';
        }
        else
        {
            result << 'N';
        }
        c++;
    }
    string result_str = result.str();
    reverse(result_str.begin(), result_str.end());
    return result_str;
}

uint32_t find_overlap_between_nodes(asg_t *g, uint32_t node1_id, uint32_t node2_id)
{
    asg_arc_t *arc_idx = asg_arc_a(g, node1_id);
    uint32_t arc_n = asg_arc_n(g, node1_id);

    for (uint32_t i = 0; i < arc_n; ++i)
    {
        uint32_t next = arc_idx[i].v;
        // cerr<<g->seq[next>>1].name<<endl;
        // cerr<<"next: "<<next<<", node2_id: "<<node2_id<<endl;
        if (next == node2_id)
        {
            return arc_idx[i].ol;
        }
    }

    return 0; 
}

void sequence_hap_chains(asg_t *g, std::vector<hap_chain_result_t> &chain_result, int *node_type)
{
    for (size_t chain_idx = 0; chain_idx < chain_result.size(); ++chain_idx)
    {
        hap_chain_result_t &result = chain_result[chain_idx];

        if (result.hap_path.empty())
        {
            result.is_valid = false;
            // cerr << "Chain " << chain_idx << ": empty path, marked invalid" << endl;
            continue;
        }
        // cerr << "Processing chain " << chain_idx << " with "
        //      << result.hap_path.size() << " nodes" << endl;
        string &seq = result.hap_sequence;
        uint32_t current_position = 0;
        for (size_t i = 0; i < result.hap_path.size(); ++i)
        {
            uint32_t node = result.hap_path[i];
            uint32_t node_id = node >> 1;
            bool is_rev = node & 1;

            string node_seq = string(g->seq[node_id].seq);
            if (is_rev)
            {
                node_seq = reverse_complement(node_seq);
            }

            result.node_positions[node] = current_position;
            if (i == 0)
            {
                seq = node_seq;
                current_position = node_seq.length();
                continue;
            }
            uint32_t prev_node = result.hap_path[i - 1];
            uint32_t prev_node_id = prev_node >> 1;
            uint32_t overlap = find_overlap_between_nodes(g, prev_node, node);

            if (overlap > 0)
            {
                if (overlap < seq.length())
                {
                    seq = seq.substr(0, seq.length() - overlap);
                    seq += node_seq;
                    current_position = seq.length();
                }
                else
                {
                    cerr << "Warning: Overlap " << overlap << " too long, direct concatenation" << endl;
                    seq += node_seq;
                    current_position = seq.length();
                }
            }
            else
            {
                const int GAP_LEN = 500;
                string gap(GAP_LEN, 'N');
                seq += gap + node_seq;
                current_position = seq.length();

                cerr << "Info: Added " << GAP_LEN << " N's between node "
                     << prev_node_id << " and " << node_id << endl;
            }
        }

        cerr << "Chain " << chain_idx << " completed: "
             << seq.length() << " bp" << endl;
    }
}

// void find_all_paths_in_hap_chains(asg_t *g, const vector<set<uint32_t>> &hap_chains, int *node_type)
// {

//     for (size_t chain_idx = 0; chain_idx < hap_chains.size(); ++chain_idx)
//     {
//         vector<hap_chain_result_t> chain_result;
//         const auto &node_set = hap_chains[chain_idx];
//         unordered_set<uint32_t> used;
//         unordered_map<uint32_t, int> indegree;
//         // if (node_set.size() < 100)
//         // {
//         //     cout << " Now < 100 hap_chains: " << chain_idx << endl;
//         //     continue;
//         // }
//         for (uint32_t u : node_set)
//         {
//             uint32_t deg = asg_arc_n(g, u);
//             asg_arc_t *arcs = asg_arc_a(g, u);
//             for (uint32_t d = 0; d < deg; ++d)
//             {
//                 if (arcs[d].del)
//                     continue;
//                 uint32_t v = arcs[d].v;
//                 if (node_set.count(v))
//                     indegree[v]++;
//             }
//         }

//         cout << "\n=== Chain " << chain_idx << " ===" << endl;
//         std::vector<uint32_t> deg_one_nodes;
//         std::vector<uint32_t> other_nodes;
//         for (uint32_t u : node_set)
//         {
//             uint32_t deg = asg_arc_n(g, u);
//             if (deg == 1)
//             {
//                 deg_one_nodes.push_back(u);
//             }
//             else
//             {
//                 other_nodes.push_back(u);
//             }
//         }

//         for (uint32_t start : deg_one_nodes)
//         {
//             if (indegree[start] > 0 || used.count(start))
//                 continue;

//             vector<uint32_t> path;
//             uint32_t cur = start;
//             used.insert(cur);
//             used.insert(start ^ 1);

//             path.push_back(cur);
//             cout << "Start: " << g->seq[start >> 1].name << endl;
//             while (true)
//             {
//                 bool extended = false;
//                 uint32_t deg = asg_arc_n(g, cur);
//                 asg_arc_t *arcs = asg_arc_a(g, cur);
//                 // if (deg == 0)
//                 // {
//                 //     uint32_t rev_deg = asg_arc_n(g, cur ^ 1);
//                 //     asg_arc_t *rev_arcs = asg_arc_a(g, cur);
//                 //     for (uint32_t d = 0; d < rev_deg; ++d)
//                 //     {
//                 //         uint32_t next = arcs[d].v;
//                 //         if (used.count(next^1))
//                 //             continue;
//                 //     }
//                 // }

//                 if (strcmp(g->seq[cur >> 1].name, "utg034464l") == 0)
//                 {
//                     // cout << "Found utg034464l and next is " << g->seq[arcs[1].v >> 1].name << endl;
//                     cur = arcs[1].v;
//                     used.insert(cur);
//                     path.push_back(cur);
//                     extended = true;
//                     continue;
//                 }
//                 for (uint32_t d = 0; d < deg; ++d)
//                 {
//                     if (arcs[d].del)
//                         continue;
//                     uint32_t next = arcs[d].v;
//                     if (deg == 2)
//                     {
//                         uint32_t next1 = arcs[0].v;
//                         uint32_t next2 = arcs[1].v;
//                         if (g->seq[next1 >> 1].len > g->seq[next2 >> 1].len)
//                         {

//                             if (node_set.count(next1) && !used.count(next1))
//                             {
//                                 // cout << "At node " << g->seq[cur >> 1].name << ", choosing longer next node " << g->seq[next1 >> 1].name << endl;
//                                 next = next1;
//                             }
//                         }
//                         else
//                         {

//                             if (node_set.count(next2) && !used.count(next2))
//                             {
//                                 // cout << "At node " << g->seq[cur >> 1].name << ", choosing longer next node " << g->seq[next2 >> 1].name << endl;
//                                 next = next2;
//                             }
//                         }
//                     }

//                     // if (deg == 2)
//                     //     next = arcs[1].v;
//                     if (!node_set.count(next))
//                         continue;
//                     if (used.count(next))
//                         continue;

//                     // 延申成功
//                     cur = next;
//                     used.insert(cur);
//                     if (asg_arc_n(g, cur) == 1 && asg_arc_n(g, cur ^ 1) == 1)
//                     {
//                         used.insert(cur ^ 1); 
//                     }
//                     // used.insert(cur ^ 1); 
//                     path.push_back(cur);
//                     cout << "Add node:  " << g->seq[cur >> 1].name << (cur & 1 ? "-" : "+") << endl;
//                     extended = true;
//                     break;
//                 }
//                 if (extended == true)
//                     continue;
//                 // TODO: add test backword
//                 cout << "start to test backword" << endl;
//                 uint32_t deg1 = asg_arc_n(g, cur ^ 1);
//                 asg_arc_t *arcs1 = asg_arc_a(g, cur ^ 1);

//                 for (uint32_t d = 0; d < deg1; ++d)
//                 {
//                     if (arcs1[d].del)
//                         continue;
//                     uint32_t next = arcs1[d].v;
//                     // if (deg == 2)
//                     //     next = arcs[1].v;
//                     if (!node_set.count(next))
//                         continue;
//                     if (used.count(next))
//                         continue;

//                     cur = next;
//                     used.insert(cur);
//                     if (asg_arc_n(g, cur) == 1 && asg_arc_n(g, cur ^ 1) == 1)
//                     {
//                         used.insert(cur ^ 1); 
//                     }
//                     // used.insert(cur ^ 1); 
//                     path.push_back(cur);
//                     cout << "Add node backward:  " << g->seq[cur >> 1].name << (cur & 1 ? "-" : "+") << endl;
//                     extended = true;
//                     break;
//                 }

//                 if (!extended)
//                     break;
//             }

//             if (path.size() > 1)
//             {
//                 cout << "Path:" << endl;
//                 for (uint32_t n : path)
//                 {
//                     cout << "  " << g->seq[n >> 1].name << (n & 1 ? "-" : "+") << endl;
//                 }

//                 hap_chain_result_t result;
//                 result.begin = path.front();
//                 result.end = path.back();
//                 result.hap_path.reserve(path.size());
//                 for (uint32_t n : path)
//                 {
//                     result.hap_path.push_back(n);
//                 }
//                 for (uint32_t n : path)
//                 {
//                     used.insert(n);
//                     used.insert(n ^ 1);
//                 }

//                 chain_result.push_back(result);
//             }
//         }

//         for (uint32_t start : other_nodes)
//         {
//             if (indegree[start] > 0 || used.count(start))
//                 continue;

//             vector<uint32_t> path;
//             uint32_t cur = start;
//             used.insert(cur);
//             used.insert(start ^ 1);

//             path.push_back(cur);
//             cout << "Start: " << g->seq[start >> 1].name << endl;
//             while (true)
//             {
//                 bool extended = false;
//                 uint32_t deg = asg_arc_n(g, cur);
//                 asg_arc_t *arcs = asg_arc_a(g, cur);
//                 // if (deg == 0)
//                 // {
//                 //     uint32_t rev_deg = asg_arc_n(g, cur ^ 1);
//                 //     asg_arc_t *rev_arcs = asg_arc_a(g, cur);
//                 //     for (uint32_t d = 0; d < rev_deg; ++d)
//                 //     {
//                 //         uint32_t next = arcs[d].v;
//                 //         if (used.count(next^1))
//                 //             continue;
//                 //     }
//                 // }

//                 // if (strcmp(g->seq[cur >> 1].name, "utg034464l") == 0)
//                 // {
//                 //     // cout << "Found utg034464l and next is " << g->seq[arcs[1].v >> 1].name << endl;
//                 //     cur = arcs[1].v;
//                 //     used.insert(cur);
//                 //     path.push_back(cur);
//                 //     extended = true;
//                 //     continue;
//                 // }
//                 for (uint32_t d = 0; d < deg; ++d)
//                 {
//                     if (arcs[d].del)
//                         continue;
//                     uint32_t next = arcs[d].v;
//                     // if (deg == 2)
//                     //     next = arcs[1].v;
//                     if (!node_set.count(next))
//                         continue;
//                     if (used.count(next))
//                         continue;

//                     // 延申成功
//                     cur = next;
//                     used.insert(cur);
//                     if (asg_arc_n(g, cur) == 1 && asg_arc_n(g, cur ^ 1) == 1)
//                     {
//                         used.insert(cur ^ 1);
//                     }
//                     // used.insert(cur ^ 1); 
//                     path.push_back(cur);
//                     cout << "Add node backward:  " << g->seq[cur >> 1].name << (cur & 1 ? "-" : "+") << endl;
//                     extended = true;
//                     break;
//                 }

//                 if (!extended)
//                     break;
//             }
//             if (path.size() > 1)
//             {
//                 cout << "Path:" << endl;
//                 for (uint32_t n : path)
//                 {
//                     cout << "  " << g->seq[n >> 1].name << (n & 1 ? "-" : "+") << endl;
//                 }

//                 hap_chain_result_t result;
//                 result.begin = path.front();
//                 result.end = path.back();
//                 result.hap_path.reserve(path.size());
//                 for (uint32_t n : path)
//                 {
//                     result.hap_path.push_back(n);
//                 }
//                 for (uint32_t n : path)
//                 {
//                     used.insert(n);
//                     used.insert(n ^ 1);
//                 }
//                 chain_result.push_back(result);
//             }
//         }

//         vector<uint32_t> unvisited_nodes;
//         for (uint32_t node : node_set)
//         {
//             if (!used.count(node))
//             {
//                 cout << "Unvisited node in chain " << chain_idx << ": " << g->seq[node >> 1].name << (node & 1 ? "-" : "+") << "  , node type:" << node_type[node] << endl;

//                 unvisited_nodes.push_back(node);
//             }
//         }
//         // TODO: add test

//         std::unordered_map<uint32_t, int> begin2chains_id;
//         std::unordered_map<uint32_t, int> end2chains_id;
//         for (int i = 0; i < chain_result.size(); i++)
//         {
//             begin2chains_id[chain_result[i].begin >> 1] = i;
//             end2chains_id[chain_result[i].end >> 1] = i;
//             // cout << "chain " << i << " " << chain_result[i].begin << " " << chain_result[i].end << " " << endl;
//         }
//         //......................................
//         std::vector<hap_chain_result_t> chain_result_2;

//         unordered_set<uint32_t> used_1;

//         unordered_set<uint32_t> used_nodes;

//         vector<vector<uint32_t>> dfs_paths;
//         unvisited_nodes_global.push_back(unvisited_nodes);
//         if (chain_result.size() > 0)
//         {
//             sequence_hap_chains(g, chain_result, node_type);
//             cerr << "Finish one chain : " << chain_result[0].hap_sequence.length() << " bp" << endl;
//             hap_results_debug.push_back(chain_result);
//         }
//     }
// }

// void find_all_paths_in_hap_chains(asg_t *g, const vector<set<uint32_t>> &hap_chains, int *node_type)
// {
//     for (size_t chain_idx = 0; chain_idx < hap_chains.size(); ++chain_idx)
//     {
//         const auto &node_set = hap_chains[chain_idx];
//         if (node_set.empty()) continue;

//         cout << "\n==============================================" << endl;
//         cout << "=== Processing Chain " << chain_idx << " (Size: " << node_set.size() << ") ===" << endl;
//         cout << "==============================================" << endl;

//         unordered_map<uint32_t, int> indegree;
//         for (uint32_t u : node_set) indegree[u] = 0; 

//         for (uint32_t u : node_set)
//         {
//             uint32_t deg = asg_arc_n(g, u);
//             asg_arc_t *arcs = asg_arc_a(g, u);
//             for (uint32_t d = 0; d < deg; ++d)
//             {
//                 if (arcs[d].del) continue;
//                 uint32_t v = arcs[d].v;
//                 if (node_set.count(v))
//                 {
//                     indegree[v]++;
//                 }
//             }
//         }

//         unordered_set<uint32_t> used;
//         vector<hap_chain_result_t> chain_result;

//         auto extract_path_forward = [&](uint32_t start) -> vector<uint32_t> {
//             vector<uint32_t> path;
//             uint32_t cur = start;

//             cout << "\n  --> [Path Start] Initiating new path from Source Node" << endl;

//             while (true)
//             {
//                 used.insert(cur);
//                 used.insert(cur ^ 1);

//                 path.push_back(cur);
//                 cout << "      [Step " << path.size() << "] Added "
//                      << "NodeID: " << cur
//                      << " | Name: " << g->seq[cur >> 1].name << (cur & 1 ? "-" : "+")
//                      << " | Len: " << g->seq[cur >> 1].len << " bp" << endl;

//                 uint32_t best_next = -1;
//                 uint32_t max_len = 0;
//                 bool found_next = false;

//                 uint32_t deg = asg_arc_n(g, cur);
//                 asg_arc_t *arcs = asg_arc_a(g, cur);

//                 for (uint32_t d = 0; d < deg; ++d)
//                 {
//                     if (arcs[d].del) continue;
//                     uint32_t next = arcs[d].v;

//                     if (!node_set.count(next) || used.count(next))
//                         continue;

//                     uint32_t next_seq_len = g->seq[next >> 1].len;
//                     if (next_seq_len > max_len)
//                     {
//                         max_len = next_seq_len;
//                         best_next = next;
//                         found_next = true;
//                     }
//                 }

//                 if (!found_next) {
//                     cout << "  --> [Path End] Reached a dead end or all branches are visited. Total nodes in this path: " << path.size() << "\n" << endl;
//                     break;
//                 }

//                 cur = best_next;
//             }
//             return path;
//         };
//         for (uint32_t start : node_set)
//         {
//             if (indegree[start] == 0 && !used.count(start) && !used.count(start ^ 1))
//             {
//                 vector<uint32_t> path = extract_path_forward(start);
//                 if (path.size() > 1) 
//                 {
//                     hap_chain_result_t result;
//                     result.begin = path.front();
//                     result.end = path.back();
//                     result.hap_path = path;
//                     chain_result.push_back(result);
//                 }
//             }
//         }

//         for (uint32_t start : node_set)
//         {
//             if (!used.count(start) && !used.count(start ^ 1))
//             {
//                 cout << "\n  [Cycle Warning] Breaking cycle. Starting from arbitrary internal node." << endl;
//                 vector<uint32_t> path = extract_path_forward(start);
//                 if (path.size() > 1)
//                 {
//                     hap_chain_result_t result;
//                     result.begin = path.front();
//                     result.end = path.back();
//                     result.hap_path = path;
//                     chain_result.push_back(result);
//                 }
//             }
//         }
//         vector<uint32_t> unvisited_nodes;
//         for (uint32_t node : node_set)
//         {
//             if (!used.count(node) && !used.count(node ^ 1))
//             {
//                 cout << "[Unvisited Warning] Leftover node in chain " << chain_idx << ": "
//                      << g->seq[node >> 1].name << (node & 1 ? "-" : "+")
//                      << " (ID: " << node << ")"
//                      << " | node_type: " << node_type[node] << endl;
//                 unvisited_nodes.push_back(node);
//             }
//         }

//         // unvisited_nodes_global.push_back(unvisited_nodes);
//         if (chain_result.size() > 0)
//         {
//             sequence_hap_chains(g, chain_result, node_type);
//             cerr << "[Success] Finished exporting one chain. Total length: " << chain_result[0].hap_sequence.length() << " bp\n" << endl;
//             // hap_results_debug.push_back(chain_result);
//         }
//     }
// }

void find_all_paths_in_hap_chains_plant(asg_t *g, const vector<set<uint32_t>> &hap_chains, int *node_type, const string &output_directory)
{
    string fa_file_path = output_directory + "/hap_chains.fa";
    string txt_file_path = output_directory + "/hap_chains_path.txt";
    ofstream fa_out(fa_file_path);
    ofstream txt_out(txt_file_path);

    if (!fa_out.is_open() || !txt_out.is_open())
    {
        cerr << "[ERROR] Cannot open output files in directory: " << output_directory << endl;
        return;
    }

    int output_counter = 1;

    for (size_t chain_idx = 0; chain_idx < hap_chains.size(); ++chain_idx)
    {
        const auto &node_set = hap_chains[chain_idx];
        if (node_set.empty())
            continue;

        cout << "\n==============================================" << endl;
        cout << "=== Processing Plant Chain " << chain_idx << " (Size: " << node_set.size() << ") ===" << endl;
        cout << "==============================================" << endl;

        unordered_map<uint32_t, int> indegree;
        unordered_map<uint32_t, int> outdegree;
        for (uint32_t u : node_set)
        {
            indegree[u] = 0;
            outdegree[u] = 0;
        }

        for (uint32_t u : node_set)
        {
            uint32_t deg = asg_arc_n(g, u);
            asg_arc_t *arcs = asg_arc_a(g, u);
            for (uint32_t d = 0; d < deg; ++d)
            {
                if (arcs[d].del)
                    continue;
                uint32_t v = arcs[d].v;
                if (node_set.count(v))
                {
                    outdegree[u]++;
                    indegree[v]++; 
                }
            }
        }

        uint32_t best_start_node = -1;
        uint32_t max_start_len = 0;
        bool found_source = false;

        for (uint32_t u : node_set)
        {
            if (indegree[u] == 0)
            {
                uint32_t len = g->seq[u >> 1].len;
                if (len > max_start_len)
                {
                    max_start_len = len;
                    best_start_node = u;
                    found_source = true;
                }
            }
        }

        if (!found_source)
        {
            cout << "  [Warning] No in-degree 0 node found (Possible cycle). Selecting the longest node overall as start." << endl;
            for (uint32_t u : node_set)
            {
                uint32_t len = g->seq[u >> 1].len;
                if (len > max_start_len)
                {
                    max_start_len = len;
                    best_start_node = u;
                }
            }
        }

        if (best_start_node == (uint32_t)-1)
            continue; 

        cout << "  --> Selected Start Node: " << best_start_node
             << " (Name: " << g->seq[best_start_node >> 1].name << (best_start_node & 1 ? "-" : "+")
             << ", Len: " << max_start_len << ")" << endl;

        vector<uint32_t> best_path;
        uint32_t global_max_path_len = 0;

        // DFS Lambda
        auto dfs_longest_path = [&](auto &self, uint32_t cur, unordered_set<uint32_t> &current_visited,
                                    vector<uint32_t> &current_path, uint32_t current_len) -> void
        {
            current_visited.insert(cur);
            current_visited.insert(cur ^ 1); 
            current_path.push_back(cur);
            current_len += g->seq[cur >> 1].len;

            bool is_leaf = true;
            uint32_t deg = asg_arc_n(g, cur);
            asg_arc_t *arcs = asg_arc_a(g, cur);

            for (uint32_t d = 0; d < deg; ++d)
            {
                if (arcs[d].del)
                    continue;
                uint32_t next = arcs[d].v;

                if (node_set.count(next) && !current_visited.count(next))
                {
                    is_leaf = false;
                    self(self, next, current_visited, current_path, current_len);
                }
            }

            if (is_leaf)
            {
                if (current_len > global_max_path_len)
                {
                    global_max_path_len = current_len;
                    best_path = current_path;
                }
            }

            current_path.pop_back();
            current_visited.erase(cur);
            current_visited.erase(cur ^ 1);
        };

        unordered_set<uint32_t> dfs_visited;
        vector<uint32_t> dfs_path;
        dfs_longest_path(dfs_longest_path, best_start_node, dfs_visited, dfs_path, 0);
        if (best_path.empty())
        {
            cout << "  [Warning] Best path is empty for chain " << chain_idx << ". Skipping..." << endl;
            continue;
        }

        unordered_set<uint32_t> final_used;
        for (uint32_t node : best_path)
        {
            final_used.insert(node);
            final_used.insert(node ^ 1);
        }

        vector<uint32_t> unvisited_nodes;
        for (uint32_t node : node_set)
        {
            if (!final_used.count(node))
            {
                cout << "[Unvisited Warning] Leftover node in chain " << chain_idx << ": "
                     << g->seq[node >> 1].name << (node & 1 ? "-" : "+")
                     << " (ID: " << node << ")"
                     << " | node_type: " << node_type[node] << endl;
                unvisited_nodes.push_back(node);
            }
        }
        unvisited_nodes_global.push_back(unvisited_nodes);

        vector<hap_chain_result_t> chain_result;
        hap_chain_result_t result;
        result.begin = best_path.front();
        result.end = best_path.back();
        result.hap_path = best_path;
        chain_result.push_back(result);
        sequence_hap_chains(g, chain_result, node_type);

        if (!chain_result.empty() && !chain_result[0].hap_path.empty() && !chain_result[0].hap_sequence.empty())
        {
            cerr << "[Success] Finished exporting plant chain. Total base pairs: " << chain_result[0].hap_sequence.length() << " bp\n"
                 << endl;

            hap_results_debug.push_back(chain_result);

            string chain_name = "hap_chains" + to_string(output_counter);

            fa_out << ">" << chain_name << "\n"
                   << chain_result[0].hap_sequence << "\n";

            txt_out << chain_name << ": ";
            for (size_t i = 0; i < best_path.size(); ++i)
            {
                uint32_t node = best_path[i];
                txt_out << g->seq[node >> 1].name << (node & 1 ? "-" : "+");
                if (i != best_path.size() - 1)
                    txt_out << " ";
            }
            txt_out << "\n";

            output_counter++; 
        }
        else
        {
            cout << "  [Warning] Generated sequence or path is empty for chain " << chain_idx << " after sequencing. Ignored." << endl;
        }
    }

    fa_out.close();
    txt_out.close();
    cout << "[INFO] All plant chains written to " << fa_file_path << " and " << txt_file_path << endl;
}

void find_all_paths_in_hap_chains(asg_t *g, const vector<set<uint32_t>> &hap_chains, int *node_type)
{
    for (size_t chain_idx = 0; chain_idx < hap_chains.size(); ++chain_idx)
    {
        const auto &node_set = hap_chains[chain_idx];
        if (node_set.empty())
            continue;

        // cout << "\n==============================================" << endl;
        // cout << "=== Processing Chain " << chain_idx << " (Size: " << node_set.size() << ") ===" << endl;
        // cout << "==============================================" << endl;

        unordered_map<uint32_t, int> indegree;
        unordered_map<uint32_t, int> outdegree;
        for (uint32_t u : node_set)
        {
            indegree[u] = 0;
            outdegree[u] = 0;
        }

        for (uint32_t u : node_set)
        {
            uint32_t deg = asg_arc_n(g, u);
            asg_arc_t *arcs = asg_arc_a(g, u);
            for (uint32_t d = 0; d < deg; ++d)
            {
                if (arcs[d].del)
                    continue;
                uint32_t v = arcs[d].v;
                if (node_set.count(v))
                {
                    outdegree[u]++; 
                    indegree[v]++;  
                }
            }
        }

        unordered_set<uint32_t> used;
        vector<hap_chain_result_t> chain_result;

        auto count_unvisited_predecessors = [&](uint32_t node) -> int
        {
            int unvisited_count = 0;
            uint32_t rev_node = node ^ 1;
            uint32_t deg = asg_arc_n(g, rev_node);
            asg_arc_t *arcs = asg_arc_a(g, rev_node);
            for (uint32_t d = 0; d < deg; ++d)
            {
                if (arcs[d].del)
                    continue;
                uint32_t rev_pred = arcs[d].v;
                uint32_t pred = rev_pred ^ 1;

                if (node_set.count(pred) && !used.count(pred))
                {
                    unvisited_count++;
                }
            }
            return unvisited_count;
        };

        auto extract_path_forward = [&](uint32_t start) -> vector<uint32_t>
        {
            vector<uint32_t> path;
            uint32_t cur = start;

            // cout << "\n  --> [Path Start] Initiating new path from Source Node" << endl;

            while (true)
            {
                used.insert(cur);
                path.push_back(cur);

                bool lock_reverse = true;

                if (indegree[cur] == 2 && outdegree[cur] == 2)
                {
                    int unv_preds = count_unvisited_predecessors(cur);
                    if (unv_preds > 0)
                    {
                        lock_reverse = false; 
                        //cout << "      [Special Branch] Node " << cur << " is a shared butterfly node. Leaving reverse path open!" << endl;
                    }
                }

                if (lock_reverse)
                {
                    used.insert(cur ^ 1);
                }

                // cout << "      [Step " << path.size() << "] Added "
                //      << "NodeID: " << cur
                //      << " | Name: " << g->seq[cur >> 1].name << (cur & 1 ? "-" : "+")
                //      << " | Len: " << g->seq[cur >> 1].len << " bp" << endl;

                uint32_t best_next = -1;
                uint32_t max_len = 0;
                bool found_next = false;

                uint32_t deg = asg_arc_n(g, cur);
                asg_arc_t *arcs = asg_arc_a(g, cur);

                for (uint32_t d = 0; d < deg; ++d)
                {
                    if (arcs[d].del)
                        continue;
                    uint32_t next = arcs[d].v;

                    if (!node_set.count(next) || used.count(next))
                        continue;

                    uint32_t next_seq_len = g->seq[next >> 1].len;
                    if (next_seq_len > max_len)
                    {
                        max_len = next_seq_len;
                        best_next = next;
                        found_next = true;
                    }
                }

                if (!found_next)
                {
                    // cout << "  --> [Path End] Reached a dead end or all branches are visited. Total nodes in this path: " << path.size() << "\n"
                    //      << endl;
                    break;
                }

                cur = best_next;
            }
            return path;
        };

        auto get_path_length = [&](const vector<uint32_t>& path) -> uint32_t
        {
            uint32_t total_len = 0;
            for (uint32_t n : path)
            {
                total_len += g->seq[n >> 1].len;
            }
            return total_len;
        };
        for (uint32_t start : node_set)
        {
            if (indegree[start] == 0 && !used.count(start) && !used.count(start ^ 1))
            {
                vector<uint32_t> path = extract_path_forward(start);
                uint32_t path_length = get_path_length(path);

                if (path_length >= 1000000) 
                {
                    hap_chain_result_t result;
                    result.begin = path.front();
                    result.end = path.back();
                    result.hap_path = path;
                    chain_result.push_back(result);
                }
            }
        }

        for (uint32_t start : node_set)
        {
            if (!used.count(start) && !used.count(start ^ 1))
            {
                vector<uint32_t> path = extract_path_forward(start);
                uint32_t path_length = get_path_length(path);

                
                if (path_length >= 1000000)
                {
                    hap_chain_result_t result;
                    result.begin = path.front();
                    result.end = path.back();
                    result.hap_path = path;
                    chain_result.push_back(result);
                }
            }
        }

        vector<uint32_t> unvisited_nodes;
        for (uint32_t node : node_set)
        {
            if (!used.count(node) && !used.count(node ^ 1))
            {
                cout << "[Unvisited Warning] Leftover node in chain " << chain_idx << ": "
                     << g->seq[node >> 1].name << (node & 1 ? "-" : "+")
                     << " (ID: " << node << ")"
                     << " | node_type: " << node_type[node] << endl;
                unvisited_nodes.push_back(node);
            }
        }
        unvisited_nodes_global.push_back(unvisited_nodes);
        if (chain_result.size() > 0)
        {
            sequence_hap_chains(g, chain_result, node_type);
            cerr << "[Success] Finished exporting one chain. Total length: " << chain_result[0].hap_sequence.length() << " bp\n"
                 << endl;
            hap_results_debug.push_back(chain_result);
        }
    }
}

std::map<size_t, std::vector<std::vector<uint32_t>>>
find_longest_paths_from_type3_nodes(
    asg_t *g,
    int *node_type,
    const std::vector<std::set<uint32_t>> &groups)
{
    std::map<size_t, std::vector<std::vector<uint32_t>>> result;

    for (size_t group_idx = 0; group_idx < groups.size(); ++group_idx)
    {
        const auto &group = groups[group_idx];
        std::vector<std::vector<uint32_t>> group_paths;

        std::cout << "\n=== Processing Group " << group_idx
                  << " (size: " << group.size() << ") ===" << std::endl;

        std::vector<uint32_t> type3_nodes;
        for (uint32_t node : group)
        {
            if (node_type[node] == 3)
            {

                uint32_t num_out = asg_arc_n(g, node);
                if (num_out != 0)
                {
                    type3_nodes.push_back(node);
                }
            }
        }

        std::cout << "Found " << type3_nodes.size()
                  << " type 3 nodes in group " << group_idx << std::endl;

        for (uint32_t start_node : type3_nodes)
        {
            std::vector<uint32_t> longest_path;
            std::set<uint32_t> permanent_visited; 
            int recursion_counter = 0;            
            const int MAX_RECURSION = 10000;      
            std::function<void(uint32_t, std::vector<uint32_t> &)>
                dfs_no_cycle = [&](uint32_t current, std::vector<uint32_t> &current_path) -> void
            {
                recursion_counter++;
                if (recursion_counter > MAX_RECURSION)
                {
                    std::cout << "WARNING: Max recursion limit (" << MAX_RECURSION
                              << ") reached for node " << g->seq[start_node >> 1].name
                              << ". Stopping DFS." << std::endl;
                    return;
                }

                current_path.push_back(current);

                if (current_path.size() > longest_path.size())
                {
                    longest_path = current_path;
                }

                uint32_t out_deg = asg_arc_n(g, current);
                asg_arc_t *out_arcs = asg_arc_a(g, current);

                for (uint32_t i = 0; i < out_deg; ++i)
                {
                    uint32_t next = out_arcs[i].v;

                    if (group.count(next) == 0 ||
                        permanent_visited.count(next) ||
                        std::find(current_path.begin(), current_path.end(), next) != current_path.end())
                        continue;
                    std::cout << "  DFS Edge: "
                              << g->seq[current >> 1].name << (current & 1 ? "-" : "+")
                              << " -> "
                              << g->seq[next >> 1].name << (next & 1 ? "-" : "+")
                              << " (depth: " << current_path.size()
                              << ", recursion: " << recursion_counter << ")" << std::endl;
                    dfs_no_cycle(next, current_path);

                    if (recursion_counter > MAX_RECURSION)
                        return;
                }

                current_path.pop_back();
            };

            std::vector<uint32_t> current_path;
            std::set<uint32_t> visited;

            if (!permanent_visited.count(start_node))
            {
                recursion_counter = 0; 
                dfs_no_cycle(start_node, current_path);

                for (uint32_t node : longest_path)
                {
                    permanent_visited.insert(node);
                }

                if (longest_path.size() > 1)
                {
                    group_paths.push_back(longest_path);

                    std::cout << "\n[Group " << group_idx << "] Longest path from type 3 node: "
                              << g->seq[start_node >> 1].name
                              << (start_node & 1 ? "-" : "+")
                              << " (length: " << longest_path.size()
                              << ", recursion calls: " << recursion_counter << ")" << std::endl;

                    std::cout << "Path: ";
                    for (size_t i = 0; i < longest_path.size(); ++i)
                    {
                        uint32_t node = longest_path[i];
                        std::cout << g->seq[node >> 1].name
                                  << (node & 1 ? "-" : "+");
                        if (i < longest_path.size() - 1)
                            std::cout << " -> ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        if (!group_paths.empty())
        {
            result[group_idx] = group_paths;
        }
    }

    return result;
}

map<uint32_t, map<uint32_t, set<uint32_t>>> *phasing_plant_version(asg_t *g, string output_directory, uint32_t **connection_forward, uint32_t **connection_backward)
{
    // define
    uint32_t **connections_count;
    uint32_t n_vtx = g->n_seq * 2;
    CALLOC(connections_count, g->n_seq);
    for (int i = 0; i < g->n_seq; i++)
    {
        CALLOC(connections_count[i], g->n_seq);
        memset(connections_count[i], 0, sizeof(*connections_count[i]));
    }
    for (int i = 0; i < g->n_seq; i++)
    {
        for (int j = 0; j < g->n_seq; j++)
        {
            connections_count[i][j] = connection_forward[i][j] + connection_backward[i][j];
            connections_count[i][j] += connection_forward[j][i] + connection_backward[j][i];
        }
    }

    int coverage[n_vtx];
    // analyze_coverage(g, coverage, n_vtx);
    cout << "start get bubbles" << endl;
    deduplicate_outgoing_arcs_safe(g);

    int node_type[n_vtx]; // Node type array
    node_type_global = (int *)calloc(n_vtx, sizeof(int));
    int node_type2[n_vtx];
    vector<bubble_t *> bubble_by_ending_begining[n_vtx];
    vector<bubble_t *> bubbles;
    Get_bubble(g, bubbles, bubble_by_ending_begining, node_type);

    std::vector<std::set<uint32_t>> hap_chains;
    std::vector<std::set<uint32_t>> filtered_components;
    auto components = find_connected_components_debug_plants(g, hap_chains, filtered_components, node_type);

    system((string("rm -r ") + output_directory).c_str());
    mkdir(output_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    node_type_global = node_type;
    hap_chains_global = hap_chains;
    find_all_paths_in_hap_chains_plant(g, hap_chains_global, node_type, output_directory);

    std::unordered_set<uint32_t> used_seq_indices;
    for (const auto &chain : hap_chains_global)
    {
        for (uint32_t node : chain)
        {
            used_seq_indices.insert(node >> 1);
        }
    }

    string scaffold_fa_path = output_directory + "/scaffold.fa";
    string scaffold_txt_path = output_directory + "/scaffold_path.txt";
    std::ofstream scaf_fa_out(scaffold_fa_path);
    std::ofstream scaf_txt_out(scaffold_txt_path);

    if (!scaf_fa_out.is_open() || !scaf_txt_out.is_open())
    {
        std::cerr << "[ERROR] Cannot open scaffold output files in: " << output_directory << std::endl;
    }
    else
    {
        int single_utg_counter = 1;

        for (uint32_t i = 0; i < (uint32_t)g->n_seq; ++i)
        {
            if (used_seq_indices.find(i) == used_seq_indices.end())
            {
                string seq_name = "scaffold_singleton_" + to_string(single_utg_counter);
                string utg_name = string(g->seq[i].name);

                string seq_str = (g->seq[i].seq) ? string(g->seq[i].seq) : "";

                scaf_fa_out << ">" << seq_name << "\n"
                            << seq_str << "\n";

                scaf_txt_out << seq_name << ": " << utg_name << "+\n";

                single_utg_counter++;
            }
        }

        scaf_fa_out.close();
        scaf_txt_out.close();
        std::cout << "[INFO] Exported " << (single_utg_counter - 1) << " unphased/single-node utgs to " << scaffold_fa_path << " and " << scaffold_txt_path << std::endl;
    }
    // =========================================================================

    exit(0);
}

std::vector<uint32_t> missing_nodes;
map<uint32_t, map<uint32_t, set<uint32_t>>> *phasing_10_7(asg_t *g, string output_directory, uint32_t **connection_forward, uint32_t **connection_backward)
{
    // define
    uint32_t **connections_count;
    uint32_t n_vtx = g->n_seq * 2;
    CALLOC(connections_count, g->n_seq);
    for (int i = 0; i < g->n_seq; i++)
    {
        CALLOC(connections_count[i], g->n_seq);
        memset(connections_count[i], 0, sizeof(*connections_count[i]));
    }
    for (int i = 0; i < g->n_seq; i++)
    {
        for (int j = 0; j < g->n_seq; j++)
        {
            connections_count[i][j] = connection_forward[i][j] + connection_backward[i][j];
            connections_count[i][j] += connection_forward[j][i] + connection_backward[j][i];
        }
    }

    int coverage[n_vtx];
    // analyze_coverage(g, coverage, n_vtx);
    cout << "start get bubbles" << endl;
    deduplicate_outgoing_arcs_safe(g);

    int node_type[n_vtx]; // Node type array
    node_type_global = (int *)calloc(n_vtx, sizeof(int));
    int node_type2[n_vtx];
    vector<bubble_t *> bubble_by_ending_begining[n_vtx];
    vector<bubble_t *> bubbles;
    Get_bubble(g, bubbles, bubble_by_ending_begining, node_type);

    std::vector<std::set<uint32_t>> hap_chains;
    std::vector<std::set<uint32_t>> filtered_components;

    // cout << "node type 4: " << endl;
    // for (uint32_t i = 0; i < n_vtx; i++)
    // {
    //     if (asg_arc_n(g, i) == 1 && asg_arc_n(g, i ^ 1) == 1)
    //     {
    //         if (node_type[i] == 0)
    //         {
    //             node_type[i] = 4;
    //             cout << g->seq[i >> 1].name << " ";
    //         }
    //     }
    // }
    // cout << endl;
    auto components = find_connected_components_debug(g, hap_chains, filtered_components, node_type);
    std::unordered_set<uint32_t> allowed_nodes;
    for (const auto &comp : components)
        for (uint32_t v : comp)
            allowed_nodes.insert(v);

    Type3Classification classification_result;
    debug_classify_type3_nodes2(g, node_type, classification_result, allowed_nodes);

    std::set<uint32_t> to_delete; 
    for (const auto &pair : classification_result.potential_middle_connections)
    {
        uint32_t a = pair.first;
        uint32_t b = pair.second;
        to_delete.insert(a);
        to_delete.insert(a ^ 1);
        to_delete.insert(b);
        to_delete.insert(b ^ 1);
    }

    std::vector<uint32_t> new_chain_middle_nodes;
    for (auto node : classification_result.chain_middle_nodes)
    {
        if (to_delete.count(node) == 0)
        {
            new_chain_middle_nodes.push_back(node);
        }
    }
    classification_result.chain_middle_nodes = std::move(new_chain_middle_nodes);

    // cout << "[Info] After cleanup, chain_middle_nodes size = " << classification_result.chain_middle_nodes.size() << endl;

    // cout << "\n==== Potential Middle Connections ====" << endl;
    // for (const auto &pair : classification_result.potential_middle_connections)
    // {
    //     uint32_t a = pair.first;
    //     uint32_t b = pair.second;
    //     cout << "Connection: " << a << " (" << g->seq[a >> 1].name << ") <--> "
    //          << b << " (" << g->seq[b >> 1].name << ")" << endl;
    // }
    // cout << "\n==== Chain Middle (Branch) Nodes ====" << endl;
    // for (auto id : classification_result.chain_middle_nodes)
    //     cout << "Node ID: " << id << ", Name: " << g->seq[id / 2].name << endl;

    // cout << "\n==== Chain End Nodes ====" << endl;
    // for (auto id : classification_result.chain_end_nodes)
    //     cout << "Node ID: " << id << ", Name: " << g->seq[id / 2].name << endl;

    add_edges_to_graph_new(g, classification_result, node_type);
    // print_graph_edges(g);
    build_bubbles2(g, bubbles_global, bubble_by_ending_begining, node_type);

    // cout << "After build_bubbles2 - bubbles_global size: " << bubbles_global.size() << endl;
    // for (int i=0; i < bubbles_global.size(); i++) {
    //     if (bubbles_global[i]) {
    //         cout << "Bubble " << i << " paths: " << bubbles_global[i]->paths.size() << endl;
    //     } else {
    //         cout << "Bubble " << i << " is NULL!" << endl;
    //     }
    // }

    hap_chains.clear();

    cout << "node type 4: " << endl;
    for (uint32_t i = 0; i < n_vtx; i++)
    {
        if (asg_arc_n(g, i) == 1 && asg_arc_n(g, i ^ 1) == 1)
        {
            if (node_type[i] == 0)
            {
                node_type[i] = 4;
                cout << g->seq[i >> 1].name << " ";
            }
        }
    }
    cout << endl;
    components = find_connected_components_debug(g, hap_chains, filtered_components, node_type);

    components1 = components;
    node_type_global = node_type;
    hap_chains_global = hap_chains;
    find_all_paths_in_hap_chains(g, hap_chains_global, node_type);
    // find_all_paths_in_hap_chains_BFS(g, hap_chains_global, node_type);
    // exit(0);
    set<uint32_t> bubble_chain_end_begin;
    map<uint32_t, uint32_t> pure_outgoing_num;
    map<uint32_t, uint32_t> pure_incoming_num;
    set<uint32_t> out_only;

    // fix_bubbles(g, bubbles, node_type2, node_type);
    analyze_bubble_node_degrees_debug(g, node_type, n_vtx, bubble_by_ending_begining,
                                      bubble_chain_end_begin, pure_outgoing_num,
                                      pure_incoming_num, out_only, components);

    extract_bubble_chains_by_dfs(g, bubble_chain_end_begin, node_type, components, bubble_chains, output_directory);

    cout << "bubble_chains.size()" << bubble_chains.size() << endl;
    // for (int i = 0; i < bubble_chains.size(); i++)
    // {
    //     cout << "Component " << i << " bubble_chains (size > 2):" << endl;
    //     for (size_t j = 0; j < bubble_chains[i].size(); ++j)
    //     {
    //         const auto &bc = bubble_chains[i][j];
    //         if (bc.nodes.size() > 2)
    //         {
    //             cout << "  Chain " << j << ", ";
    //             cout << g->seq[bc.begin >> 1].name << " -> " << g->seq[bc.end >> 1].name << " : ";
    //             for (auto node : bc.nodes)
    //             {
    //                 cout << g->seq[node >> 1].name << " ";
    //             }
    //             cout << endl;
    //         }
    //     }
    // }
    std::map<uint32_t, std::map<uint32_t, std::set<uint32_t>>> *bubble_chain_begin_end_nodes_buf =
        new std::map<uint32_t, std::map<uint32_t, std::set<uint32_t>>>;
    collect_bubble_chains_to_map(bubble_chains, bubble_chain_begin_end_nodes_buf);

    // cout << "==== bubble_chain_begin_end_nodes_buf Graph ====" << endl;
    // int bubble_chain_begin_end_nodes_buf_size = 0;
    // for (auto &begin_end_pair : *bubble_chain_begin_end_nodes_buf)
    // {
    //     uint32_t bubble_chain_begin = begin_end_pair.first;
    //     cout << bubble_chain_begin << " Bubble Chain Begin/End Node: " << g->seq[bubble_chain_begin / 2].name << endl;

    //     for (auto &target_map : begin_end_pair.second)
    //     {
    //         uint32_t target_node = target_map.first; 
    //         cout << target_node << "  Connected Target Node: " << g->seq[target_node / 2].name << endl;

    //         size_t set_size = target_map.second.size();
    //         cout << "    Number of Connected Nodes from " << g->seq[target_node / 2].name << ": " << set_size << endl;
    //         bubble_chain_begin_end_nodes_buf_size++;
    //         cout << "    Connected Nodes from " << g->seq[target_node / 2].name << ": ";
    //         for (uint32_t connected_node : target_map.second)
    //         {
    //             cout << g->seq[connected_node >> 1].name << " ";
    //         }
    //         cout << endl;
    //     }
    // }
    // cout << "bubble_chain_begin_end_nodes_buf_size:  " << bubble_chain_begin_end_nodes_buf_size << endl;
    // cout << "==== End of bubble_chain_begin_end_nodes_buf Graph ====" << endl;

    for (auto a : *bubble_chain_begin_end_nodes_buf)
    {
        for (auto b : a.second)
        {
            for (auto bub : bubble_by_ending_begining[a.first])
            {
                if (bub->begNode == a.first)
                {
                    for (auto node : bub->starting_arcs)
                    {
                        (*bubble_chain_begin_end_nodes_buf)[a.first][b.first].insert(node);
                    }
                }
            }
            for (auto bub : bubble_by_ending_begining[b.first])
            {
                if (bub->endNode == b.first)
                {
                    for (auto node : bub->ending_arcs)
                    {
                        (*bubble_chain_begin_end_nodes_buf)[a.first][b.first].insert(node);
                    }
                }
            }
        }
    }
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes_buf2 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();
    for (auto a : *bubble_chain_begin_end_nodes_buf)
    {
        for (auto b : a.second)
        {
            if ((*bubble_chain_begin_end_nodes_buf2).find(a.first) == (*bubble_chain_begin_end_nodes_buf2).end())
            {
                (*bubble_chain_begin_end_nodes_buf2)[a.first] = map<uint32_t, set<uint32_t>>();
            }
            if ((*bubble_chain_begin_end_nodes_buf2).find(b.first ^ 1) == (*bubble_chain_begin_end_nodes_buf2).end())
            {
                (*bubble_chain_begin_end_nodes_buf2)[b.first ^ 1] = map<uint32_t, set<uint32_t>>();
            }
            (*bubble_chain_begin_end_nodes_buf2)[a.first][b.first] = b.second;
            (*bubble_chain_begin_end_nodes_buf2)[b.first ^ 1][a.first ^ 1] = b.second;
        }
    }
    delete bubble_chain_begin_end_nodes_buf;

    int self_loop_count = 0;
    vector<string> self_loop_nodes;
    for (auto &[begin_node, end_map] : *bubble_chain_begin_end_nodes_buf2)
    {
        for (auto &[end_node, _] : end_map)
        {
            if ((begin_node >> 1) == (end_node >> 1))
            {
                self_loop_count++;

                string node_name = g->seq[begin_node >> 1].name;
                self_loop_nodes.push_back(node_name);
            }
        }
    }

    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes_buf3 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_end_begin_nodes_buf3 = new map<uint32_t, map<uint32_t, set<uint32_t>>>();
    for (auto a : *bubble_chain_begin_end_nodes_buf2)
    {
        // if (unaccessable.find(a.first >> 1) == unaccessable.end())
        {
            for (auto b : a.second)
            {
                // if (unaccessable.find(b.first >> 1) == unaccessable.end() && b.second.size() > 1)
                {
                    if ((*bubble_chain_begin_end_nodes_buf3).find(a.first) == (*bubble_chain_begin_end_nodes_buf3).end())
                    {
                        (*bubble_chain_begin_end_nodes_buf3)[a.first] = map<uint32_t, set<uint32_t>>();
                    }
                    if ((*bubble_chain_begin_end_nodes_buf3).find(b.first ^ 1) == (*bubble_chain_begin_end_nodes_buf3).end())
                    {
                        (*bubble_chain_begin_end_nodes_buf3)[b.first ^ 1] = map<uint32_t, set<uint32_t>>();
                    }
                    (*bubble_chain_begin_end_nodes_buf3)[a.first][b.first] = b.second;
                    (*bubble_chain_begin_end_nodes_buf3)[b.first ^ 1][a.first ^ 1] = b.second;
                    if ((*bubble_chain_end_begin_nodes_buf3).find(b.first) == (*bubble_chain_end_begin_nodes_buf3).end())
                    {
                        (*bubble_chain_end_begin_nodes_buf3)[b.first] = map<uint32_t, set<uint32_t>>();
                    }
                    if ((*bubble_chain_end_begin_nodes_buf3).find(a.first ^ 1) == (*bubble_chain_end_begin_nodes_buf3).end())
                    {
                        (*bubble_chain_end_begin_nodes_buf3)[a.first ^ 1] = map<uint32_t, set<uint32_t>>();
                    }
                    (*bubble_chain_end_begin_nodes_buf3)[b.first][a.first] = b.second;
                    (*bubble_chain_end_begin_nodes_buf3)[a.first ^ 1][b.first ^ 1] = b.second;
                }
            }
        }
    }
    delete bubble_chain_begin_end_nodes_buf2;
    cout << "Number of self-loop nodes: " << self_loop_count << endl;
    cout << "Self-loop node names:" << endl;
    for (const auto &name : self_loop_nodes)
    {
        cout << name << endl;
    }
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_begin_end_nodes = new map<uint32_t, map<uint32_t, set<uint32_t>>>();

    for (const auto &entry : *bubble_chain_begin_end_nodes_buf3)
    {
        uint32_t begin_node = entry.first;
        const auto &end_map = entry.second;

        for (const auto &end_entry : end_map)
        {
            uint32_t end_node = end_entry.first;
            const auto &node_set = end_entry.second;

            uint64_t all_len = 0;
            for (const auto &n : node_set)
            {
                all_len += g->seq[n >> 1].len;
            }

            if (all_len > 50000)
            {
                (*bubble_chain_begin_end_nodes)[begin_node][end_node] = node_set;

                // begin_end_to_group[{begin_node, end_node}] = group_idx;
            }
        }
    }

    for (auto &begin_end_pair : *bubble_chain_begin_end_nodes)
    {
        for (auto &end_pair : begin_end_pair.second)
        {
            std::set<uint32_t> new_set;

            for (uint32_t node : end_pair.second)
            {
                if ((begin_end_pair.first >> 1) == (node >> 1) ||
                    (end_pair.first >> 1) == (node >> 1))
                {
                    new_set.insert(node >> 1); 
                }
                else
                {
                    new_set.insert(node >> 1);
                }
            }

            end_pair.second = std::move(new_set);
        }
    }

    uint32_t branch_counter = 0;
    std::set<uint32_t> all_nodes;
    std::set<uint32_t> bubble_nodes;
    std::set<uint32_t> hap_nodes;
    std::set<uint32_t> bubble_chain_nodes;


    for (auto &a : *bubble_chain_begin_end_nodes)
    {
        for (auto &b : a.second)
        {
            branch_counter++;
            bubble_nodes.insert(b.second.begin(), b.second.end());
            all_nodes.insert(b.second.begin(), b.second.end());
        }
    }

    for (auto &a : hap_chains)
    {
        for (auto &n : a)
        {
            hap_nodes.insert(n >> 1);
            all_nodes.insert(n >> 1);
        }
    }

    uint32_t bubble_chain_count = 0;

    for (auto &chain_group : bubble_chains)    {
        for (auto &bc : chain_group)
        {
            bubble_chain_count++;

            for (auto &n : bc.nodes)
            {
                bubble_chain_nodes.insert(n >> 1);
                all_nodes.insert(n >> 1);
            }
        }
    }

    // find_longest_paths_from_type3_nodes(g, node_type, components1);
    std::cout << "Nodes in original Graph: " << g->n_seq << std::endl;
    std::cout << "Total Nodes in Chain Graph: " << all_nodes.size() << std::endl;
    std::cout << "Branches: " << branch_counter << std::endl;
    std::cout << "Bubble Nodes: " << bubble_nodes.size() << std::endl;
    std::cout << "Hap Nodes: " << hap_nodes.size() << std::endl;
    std::cout << "Bubble Chains: " << bubble_chain_count << std::endl;
    std::cout << "Unique Nodes in Bubble Chains: " << bubble_chain_nodes.size() << std::endl;


    for (uint32_t i = 0; i < g->n_seq; i++)
    {
        if (all_nodes.find(i) == all_nodes.end())
        {
            if (std::find(missing_nodes.begin(), missing_nodes.end(), i >> 1) == missing_nodes.end())
            {
                missing_nodes.push_back(i >> 1);
            }
        }
    }
    uint32_t missing_length = 0;
    std::cout << "Missing Nodes (" << missing_nodes.size() << "): ";
    for (size_t i = 0; i < missing_nodes.size(); i++)
    {
        std::cout << g->seq[missing_nodes[i] >> 1].name;
        missing_length += g->seq[missing_nodes[i] >> 1].len;
        if (i + 1 < missing_nodes.size())
            std::cout << ",";
    }
    std::cout << std::endl;
    std::cout << "Total Missing Length: " << missing_length << " bp" << std::endl;
    return bubble_chain_begin_end_nodes;
}

#endif
