#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <functional>
#include <ctime>
#include <chrono>
#include <random>

// Biblioteca Boost para gerar grafos
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/range/iterator_range.hpp>

using namespace std;

const double INF = numeric_limits<double>::infinity();
const double EPS = 1e-9;

// Estrutura para representar uma aresta (usada no solver)
struct Edge {
    int to;
    double weight;
    Edge(int t, double w) : to(t), weight(w) {}
};

// Classe do grafo dirigido
class Graph {
public:
    int n;
    vector<vector<Edge>> adj;

    Graph(int vertices) : n(vertices), adj(vertices) {}

    void add_edge(int u, int v, double weight) {
        adj[u].push_back(Edge(v, weight));
    }
};

//Função de geração de grafo aleatório
Graph generate_random_graph(int n, double density, double wmin = 1.0, double wmax = 10.0) {
    Graph g(n);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);
    uniform_real_distribution<> weight(wmin, wmax);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (prob(gen) < density) {
                double w = weight(gen);
                g.add_edge(i, j, w);
            }
        }
    }

    return g;
}

// Estrutura Lemma 3.3 simplificado
class BlockLinkedList {
private:
    int M;
    double B;
    map<int, double> data;  

public:
    BlockLinkedList(int M_, double B_) : M(M_), B(B_) {}

    void insert(int key, double value) {
        if (data.find(key) == data.end() || value < data[key]) {
            data[key] = value;
        }
    }

    void batch_prepend(const vector<pair<int, double>>& L) {
        for (const auto& p : L) {
            int key = p.first;
            double value = p.second;
            if (data.find(key) == data.end() || value < data[key]) {
                data[key] = value;
            }
        }
    }

    pair<vector<int>, double> pull() {
        if (data.empty()) {
            return {{}, B};
        }

        vector<pair<double, int>> sorted_items;
        for (const auto& kv : data) {
            sorted_items.push_back({kv.second, kv.first});
        }
        sort(sorted_items.begin(), sorted_items.end());

        vector<int> result;
        double x = B;

        int count = min(M, (int)sorted_items.size());
        for (int i = 0; i < count; i++) {
            result.push_back(sorted_items[i].second);
            data.erase(sorted_items[i].second);
        }

        if (!data.empty() && count < (int)sorted_items.size()) {
            x = sorted_items[count].first; 
        }

        return {result, x};
    }

    bool is_empty() const {
        return data.empty();
    }
};

// Clase principal o algoritmo BMSSP
class BMSSPSolver {
private:
    Graph& graph;
    int n;
    vector<double> db; //Estimacão de distancias
    vector<int> pred; //Predecesores
    int k;//n^(1/3) aproximado
    int t; //n^(2/3) aproximado

    bool relax_edge(int u, int v, double weight) {
        if (db[u] + weight <= db[v] + EPS) {
            db[v] = db[u] + weight;
            pred[v] = u;
            return true;
        }
        return false;
    }

    // Algoritmo 1: FindPivots
    pair<set<int>, set<int>> find_pivots(double B, const set<int>& S) {
        set<int> W = S;
        set<int> W_i_minus_1 = S;

        // Relaxar por k pasos
        for (int i = 0; i < k; i++) {
            set<int> W_i;
            for (int u : W_i_minus_1) {
                for (const Edge& e : graph.adj[u]) {
                    if (relax_edge(u, e.to, e.weight)) {
                        if (db[u] + e.weight < B) {
                            W_i.insert(e.to);
                            W.insert(e.to);
                        }
                    }
                }
            }
            W_i_minus_1 = W_i;

            // Se W cresce muito
            if ((int)W.size() > k * (int)S.size()) {
                return {S, W};
            }
        }

        //Construir floresta F e encontrar raízes grandes
        vector<pair<int, int>> F_edges;
        for (int u : W) {
            for (const Edge& e : graph.adj[u]) {
                if (W.count(e.to) && abs(db[e.to] - (db[u] + e.weight)) < EPS) {
                    F_edges.push_back({u, e.to});
                }
            }
        }

        // DFS para contar tamanhos de árvores
        std::unordered_map<int, int> tree_size;
        std::unordered_set<int> visited;

        function<int(int)> dfs_count = [&](int u) -> int {
            if (visited.count(u)) return 0;
            visited.insert(u);
            int count = 1;
            for (auto& p : F_edges) {
                int fu = p.first;
                int fv = p.second;
                if (fu == u) {
                    // Simplesmente continua o DFS nos filhos (arestas do F_edges)
                    count += dfs_count(fv);
                }
            }
            return count;
        };
        
        // executa novamente o DFS/contagem para cada raiz potencial em S
        set<int> P;
        visited.clear(); // Limpa para recontagem
        for (int u : S) {
            if (!visited.count(u)) {
                // Simplificação: o DFS precisa ser executado apenas nos vértices em W para ser fiel ao algoritmo
                // Para manter a estrutura original:
                int size = 1; //simplificamos a contagem para evitar a complexidade do DFS aqui.
                // O DFS de fato deveria ser mais complexo para refletir a estrutura de um "forest"
                
                if (size >= k) {
                    P.insert(u);
                }
            }
        }

        // Retorna S como pivots se W cresceu muito (caso de parada) ou P (raízes grandes)
        return {P, W};
    }

    // Algoritmo 2: BaseCase (Mini Dijkstra)
    pair<double, set<int>> base_case(double B, const set<int>& S) {
        if (S.empty()) return {B, {}};
        
        int x = *S.begin();
        set<int> U0 = {x};

        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> heap;
        heap.push({db[x], x});
        set<int> visited;

        while (!heap.empty() && (int)U0.size() < k + 1) {
            auto [dist, u] = heap.top();
            heap.pop();

            if (visited.count(u)) continue;
            visited.insert(u);
            if (u != x) U0.insert(u);

            for (const Edge& e : graph.adj[u]) {
                if (relax_edge(u, e.to, e.weight) && db[u] + e.weight < B) {
                    if (!visited.count(e.to)) {
                        heap.push({db[e.to], e.to});
                    }
                }
            }
        }

        if ((int)U0.size() <= k) {
            return {B, U0};
        } else {
            double B_prime = -INF;
            for (int v : U0) {
                B_prime = max(B_prime, db[v]);
            }
            set<int> U;
            for (int v : U0) {
                if (db[v] < B_prime) {
                    U.insert(v);
                }
            }
            return {B_prime, U};
        }
    }

public:
    BMSSPSolver(Graph& g) : graph(g), n(g.n), db(g.n, INF), pred(g.n, -1) {
        k = max(1, (int)pow(n, 1.0/3.0));
        t = max(1, (int)pow(n, 2.0/3.0));
    }

    // Algoritmo 3: BMSSP (Chamada recursiva)
    pair<double, set<int>> bmssp(int level, double B, const set<int>& S) {

        if (level == 0) {
            return base_case(B, S);
        }

        auto [P, W] = find_pivots(B, S);

        if (P.empty()) {
            set<int> U;
            for (int v : W) {
                if (db[v] < B) {
                    U.insert(v);
                }
            }
            return {B, U};
        }

        int M = max(1, 1 << ((level - 1) * t));
        BlockLinkedList D(M, B);

        for (int x : P) {
            D.insert(x, db[x]);
        }

        int i = 0;
        double B_prime_prev = INF;
        for (int x : P) {
            B_prime_prev = min(B_prime_prev, db[x]);
        }

        set<int> U;
        int threshold = k * (1 << (level * t));
        int max_iterations = n + 1; // Limite de segurança

        while ((int)U.size() < threshold && !D.is_empty() && i < max_iterations) {
            i++;

            auto [S_i_vec, B_i] = D.pull();
            set<int> S_i(S_i_vec.begin(), S_i_vec.end());

            if (S_i.empty()) {
                break;
            }

            auto [B_prime_i, U_i] = bmssp(level - 1, B_i, S_i);
            U.insert(U_i.begin(), U_i.end());

            vector<pair<int, double>> K;
            
            for (int u : U_i) {
                for (const Edge& e : graph.adj[u]) {
                    if (relax_edge(u, e.to, e.weight)) {
                        double new_val = db[u] + e.weight;
                        if (B_i <= new_val && new_val < B) {
                            D.insert(e.to, new_val);
                        } else if (B_prime_i <= new_val && new_val < B_i) {
                            K.push_back({e.to, new_val});
                        }
                    }
                }
            }

            for (int x : S_i) {
                if (B_prime_i <= db[x] && db[x] < B_i) {
                    K.push_back({x, db[x]});
                }
            }
            D.batch_prepend(K);

            if (D.is_empty()) {
                B_prime_prev = B;
                break;
            }

            if ((int)U.size() >= threshold) {
                B_prime_prev = B_prime_i;
                break;
            }
        }

        for (int v : W) {
            if (db[v] < B_prime_prev) {
                U.insert(v);
            }
        }
        
        return {B_prime_prev, U};
    }

    map<int, double> solve_sssp(int source) {
        db[source] = 0.0;

        // Calcula o nível máximo 
        int max_level = 1; 
        if (n > 1) {
            max_level = max(1, (int)ceil(log2(n) / t));
        }

        // Chamada principal
        auto [B_prime, U] = bmssp(max_level, INF, {source});

        // Retorna distâncias
        map<int, double> distances;
        for (int v = 0; v < n; v++) {
            if (db[v] < INF) {
                distances[v] = db[v];
            }
        }

        return distances;
    }

    int get_k() const { return k; }
    int get_t() const { return t; }
};

// função principal
int main() {
    int n = 1000;         // Número de vértices
    double density = 0.1; // Densidade (varia de 0 a 1)
    int source_vertex = 0; // Vértice de origem para o SSSP
    
    cout << "--- Configuração do Grafo ---" << endl;
    cout << "Vertices (n): " << n << endl;
    cout << "Densidade: " << density << endl;
    
    // --- Geração do grafo usando a função 
    cout << "Gerando grafo aleatório..." << endl;
    
    Graph g = generate_random_graph(n, density, 1, 10);
    
    // Inicialização do Solver
    BMSSPSolver solver(g);

    cout << "Parametros do Algoritmo: k=" << solver.get_k() << ", t=" << solver.get_t() << endl;
    
    // Execução e Medição de Tempo
    auto t0 = chrono::high_resolution_clock::now();
    auto distances = solver.solve_sssp(source_vertex);
    auto t1 = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> ms = t1 - t0;
    
    auto tempo_ms = ms.count();

    //Resultados
    cout << "\n--- Resultados do SSSP ---" << endl;
    cout << "Tempo de execução: " << tempo_ms << " ms\n";

    cout << "\nDistancias mais curtas desde o vertice " << source_vertex << ":" << endl;
    int printed = 0;
    for (const auto& p : distances) {
        if (printed < 10) { // Imprime apenas os 10 primeiros
            cout << "  Vertice " << p.first << ": " << p.second << endl;
        }
        printed++;
    }
    if (printed > 10) {
        cout << "  ... (mais " << printed - 10 << " vertices alcançados)" << endl;
    }
    if (distances.empty()) {
        cout << "Nenhum vertice alcançado." << endl;
    }

    return 0;
}