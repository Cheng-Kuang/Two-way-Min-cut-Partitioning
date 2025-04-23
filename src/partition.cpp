#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <chrono>

using namespace std;

const double MAX_double = 1.79769e+308;
chrono::steady_clock::time_point gStartTime;
const long TIMEOUT_THRESHOLD = 150; 


bool timeoutReached() {
    auto now = chrono::steady_clock::now();
    return chrono::duration_cast<chrono::seconds>(now - gStartTime).count() >= TIMEOUT_THRESHOLD;
}

struct Cell {
    string name;
    long long area;
};

struct Die {
    long long remain_area;
};

struct connectedCell {
    string name;
    double weight;
};

struct realNet {
    unordered_map<string, unordered_map<string, vector<string>>> connected_by_Nets;
    unordered_map<string, long long> Net_realWeight;
    unordered_map<string, bool> Net_beenUsed;
};


struct CellInfo {
    string name;
    long long areaA;
    long long areaB;
};


struct PartitionResult {
    unordered_map<string, bool> indieA_map;
    vector<string> cellsA;
    vector<string> cellsB;
    Die dieA;
    Die dieB;
    double cost;
};
PartitionResult gPartNet; 

void readfile(const string& filename,
    pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    Die& dieA, Die& dieB,
    unordered_map<string, vector<connectedCell>>& nets,
    realNet& realnet,
    vector<CellInfo>& cellInfos) {

    ifstream infile(filename);
    string line;
    long long max_cell_areaA = 0;
    long long max_cell_areaB = 0;

    if (!infile) {
        cerr << "Error: Unable to open file!" << endl;
        return;
    }

    long long dieWidth = 0, dieHeight = 0;

    while (getline(infile, line)) {
        stringstream ss(line);
        string keyword;
        ss >> keyword;

        if (keyword == "NumTechs") {
            continue;
        }
        else if (keyword == "Tech") {
            string techName;
            long long numLibCells;
            ss >> techName >> numLibCells;
            for (long long i = 0; i < numLibCells; i++) {
                getline(infile, line);
                stringstream libss(line);
                string libKeyword, libName;
                long long width, height;
                libss >> libKeyword >> libName >> width >> height;
                if (libKeyword == "LibCell") {
                    if (techName == "TA") {
                        libCells.first[libName] = width * height;
                        max_cell_areaA = max(width * height, max_cell_areaA);
                    }
                    else {
                        libCells.second[libName] = width * height;
                        max_cell_areaB = max(width * height, max_cell_areaB);
                    }
                }
            }
        }
        else if (keyword == "DieSize") {
            ss >> dieWidth >> dieHeight;
        }
        else if (keyword == "DieA" || keyword == "DieB") {
            int utilization;
            string tech;
            ss >> tech >> utilization;
            if (keyword == "DieA") {
                dieA.remain_area = dieWidth * dieHeight * utilization / 100;
            }
            else {
                dieB.remain_area = dieWidth * dieHeight * utilization / 100;
            }
        }
        else if (keyword == "NumCells") {
            long long numCells;
            ss >> numCells;
            for (long long i = 0; i < numCells; ++i) {
                getline(infile, line);
                stringstream cellss(line);
                string cellKeyword, cellName, libName;
                cellss >> cellKeyword >> cellName >> libName;

                CellInfo ci;
                ci.name = cellName;
                ci.areaA = libCells.first[libName];
                ci.areaB = libCells.second[libName];
                cellInfos.push_back(ci);

                libCells.first[cellName] = libCells.first[libName];
                libCells.second[cellName] = libCells.second[libName];
            }
        }
        else if (keyword == "NumNets") {
            long long numNets;
            unordered_map<string, unordered_map<string, double>> netConnections;
            ss >> numNets;
            for (long long i = 0; i < numNets; i++) {
                getline(infile, line);
                stringstream netss(line);
                string netKeyword, netName;
                long long numPins;
                double weight;
                netss >> netKeyword >> netName >> numPins >> weight;
                realnet.Net_realWeight[netName] = weight;   // record real weight for cutsize calculate
                realnet.Net_beenUsed[netName] = false;
                if (netKeyword == "Net") {
                    vector<string> connectedCells;
                    for (long long j = 0; j < numPins; j++) {
                        getline(infile, line);
                        stringstream cellss(line);
                        string cellKeyword, cellName;
                        cellss >> cellKeyword >> cellName;
                        if (cellKeyword == "Cell") {
                            connectedCells.push_back(cellName);
                        }
                    }
                    double connectionWeight = weight / (numPins - 1);
                    for (long long j = 0; j < numPins; j++) {
                        for (long long k = j + 1; k < numPins; k++) {
                            netConnections[connectedCells[j]][connectedCells[k]] += connectionWeight;
                            netConnections[connectedCells[k]][connectedCells[j]] += connectionWeight;
                            realnet.connected_by_Nets[connectedCells[j]][connectedCells[k]].push_back(netName);
                            realnet.connected_by_Nets[connectedCells[k]][connectedCells[j]].push_back(netName);
                        }
                    }
                }
            }
            for (auto it = netConnections.begin(); it != netConnections.end(); it++) {
                for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                    nets[it->first].push_back({ it2->first, it2->second });
                }
            }
        }
    }
}


double computeInitialCost(const unordered_map<string, bool>& indieA_map, const unordered_map<string, vector<connectedCell>>& nets) {
    double cost = 0;
    for (const auto& entry : indieA_map) {
        if (entry.second) {
            const string& cellA = entry.first;
            if (nets.find(cellA) != nets.end()) {
                for (const auto& conn : nets.at(cellA)) {
                    if (indieA_map.find(conn.name) != indieA_map.end() && !indieA_map.at(conn.name)) {
                        cost += conn.weight;
                    }
                }
            }
        }
    }
    return cost;
}

PartitionResult correct_the_partition(PartitionResult result, const vector<CellInfo>& cellInfos,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {
    if (result.dieA.remain_area >= 0 && result.dieB.remain_area >= 0)
        return result;

    unordered_map<string, CellInfo> cellMap;
    for (const auto& cell : cellInfos)
        cellMap[cell.name] = cell;


    while (result.dieA.remain_area < 0 || result.dieB.remain_area < 0) {
        vector<pair<string, long long>> candB;
        vector<pair<string, long long>> candA;
        for (const auto& entry : result.indieA_map) {
            const CellInfo& ci = cellMap[entry.first];
            if (!entry.second) {
                long long improve = ci.areaB - ci.areaA;
                candB.push_back({ entry.first, improve });
            }
            else {
                long long improve = ci.areaA - ci.areaB;
                candA.push_back({ entry.first, improve });
            }
        }

        sort(candB.begin(), candB.end(), [](const pair<string, long long>& a, const pair<string, long long>& b) {
            return a.second > b.second;
            });
        sort(candA.begin(), candA.end(), [](const pair<string, long long>& a, const pair<string, long long>& b) {
            return a.second > b.second;
            });

        // SWAP
        long long idxA = 0, idxB = 0;
        while ((result.dieA.remain_area < 0 || result.dieB.remain_area < 0) && idxA < candA.size() && idxB < candB.size()) {
            if(timeoutReached())
                return gPartNet;

            if (result.dieA.remain_area < 0 && result.dieB.remain_area < 0) {
                const CellInfo& ciB = cellMap[candB[idxB].first];
                result.indieA_map[candB[idxB].first] = true;
                result.dieA.remain_area -= ciB.areaA;
                result.dieB.remain_area += ciB.areaB;

                const CellInfo& ciA = cellMap[candA[idxA].first];
                result.indieA_map[candA[idxA].first] = false;
                result.dieA.remain_area += ciA.areaA;
                result.dieB.remain_area -= ciA.areaB;
                idxA++;
                idxB++;
            }
            else if (result.dieA.remain_area < 0) {
                const CellInfo& ciA = cellMap[candA[idxA].first];
                result.indieA_map[candA[idxA].first] = false;
                result.dieA.remain_area += ciA.areaA;
                result.dieB.remain_area -= ciA.areaB;
                idxA++;
            }
            else if (result.dieB.remain_area < 0) {
                const CellInfo& ciB = cellMap[candB[idxB].first];
                result.indieA_map[candB[idxB].first] = true;
                result.dieA.remain_area -= ciB.areaA;
                result.dieB.remain_area += ciB.areaB;
                idxB++;
            }
        }
    }

    result.cellsA.clear();
    result.cellsB.clear();
    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result.cost = (result.dieA.remain_area >= 0 && result.dieB.remain_area >= 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;

    return result;
}

PartitionResult initialPartition_mustok(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;


    for (const auto& cell : cellInfos) {
        result.indieA_map[cell.name] = false;
        result.dieB.remain_area -= cell.areaB;
    }


    vector<pair<string, long long>> candidate;
    for (const auto& cell : cellInfos) {
        long long diff = cell.areaB - cell.areaA;
        candidate.push_back({ cell.name, diff });
    }
    sort(candidate.begin(), candidate.end(), [](const pair<string, long long>& a, const pair<string, long long>& b) {
        return a.second > b.second;
        });


    unordered_map<string, CellInfo> cellMap;
    for (const auto& cell : cellInfos)
        cellMap[cell.name] = cell;

    for (const auto& entry : candidate) {
        if (result.dieB.remain_area >= 0)
            break;
        if (entry.second > 0) {
            const CellInfo& ci = cellMap[entry.first];
            if (result.dieA.remain_area >= ci.areaA) {
                result.indieA_map[ci.name] = true;
                result.dieA.remain_area -= ci.areaA;
                result.dieB.remain_area += ci.areaB;
            }
        }
    }

    for (const auto& cell : cellInfos) {
        if(timeoutReached())
            return gPartNet;

        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }
    result.cost = (result.dieA.remain_area > 0 && result.dieB.remain_area > 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_sortA(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;

    vector<CellInfo> sortedCells = cellInfos;
    sort(sortedCells.begin(), sortedCells.end(), [](const CellInfo& a, const CellInfo& b) {
        return a.areaA < b.areaA;
        });

    for (const auto& cell : sortedCells) {
        if (result.dieA.remain_area >= sortedCells[sortedCells.size() - 1].areaA) {
            result.indieA_map[cell.name] = true;
            result.dieA.remain_area -= cell.areaA;
        }
        else {
            result.indieA_map[cell.name] = false;
            result.dieB.remain_area -= cell.areaB;
        }
    }

    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result = correct_the_partition(result, cellInfos, libCells, nets);

    result.cost = (result.dieA.remain_area >= 0 && result.dieB.remain_area >= 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_random(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;

    vector<CellInfo> randomCells = cellInfos;
    random_shuffle(randomCells.begin(), randomCells.end());

    for (const auto& cell : randomCells) {
        double coin = static_cast<double>(rand()) / RAND_MAX;
        bool assigned = false;
        if (coin < 0.5) {
            if (result.dieA.remain_area >= cell.areaA) {
                result.indieA_map[cell.name] = true;
                result.dieA.remain_area -= cell.areaA;
                assigned = true;
            }
            else if (result.dieB.remain_area >= cell.areaB) {
                result.indieA_map[cell.name] = false;
                result.dieB.remain_area -= cell.areaB;
                assigned = true;
            }
        }
        else {
            if (result.dieB.remain_area >= cell.areaB) {
                result.indieA_map[cell.name] = false;
                result.dieB.remain_area -= cell.areaB;
                assigned = true;
            }
            else if (result.dieA.remain_area >= cell.areaA) {
                result.indieA_map[cell.name] = true;
                result.dieA.remain_area -= cell.areaA;
                assigned = true;
            }
        }
        if (!assigned) {
            if (result.dieA.remain_area >= cell.areaA) {
                result.indieA_map[cell.name] = true;
                result.dieA.remain_area -= cell.areaA;
            }
            else {
                result.indieA_map[cell.name] = false;
                result.dieB.remain_area -= cell.areaB;
            }
        }
    }

    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    if (result.dieA.remain_area < 0 || result.dieB.remain_area < 0)
        result = correct_the_partition(result, cellInfos, libCells, nets);

    result.cost = (result.dieA.remain_area >= 0 && result.dieB.remain_area >= 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_greedy(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;
    for (const auto& cell : cellInfos) {
        if (result.dieA.remain_area >= cell.areaA && cell.areaA < cell.areaB) {
            result.indieA_map[cell.name] = true;
            result.dieA.remain_area -= cell.areaA;
        }
        else {
            result.indieA_map[cell.name] = false;
            result.dieB.remain_area -= cell.areaB;
        }
    }
    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result = correct_the_partition(result, cellInfos, libCells, nets);

    result.cost = (result.dieA.remain_area > 0 && result.dieB.remain_area > 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_original(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;
    long long maxA = 0;
    for (const auto& cell : cellInfos) {
		maxA = max(maxA, cell.areaA);
        if (result.dieA.remain_area >= 3*maxA) {
            result.indieA_map[cell.name] = true;
            result.dieA.remain_area -= cell.areaA;
        }
        else {
            result.indieA_map[cell.name] = false;
            result.dieB.remain_area -= cell.areaB;
        }
    }
    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result = correct_the_partition(result, cellInfos, libCells, nets);

    result.cost = (result.dieA.remain_area > 0 && result.dieB.remain_area > 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_73(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;

    for (const auto& cell : cellInfos) {
        if (result.dieA.remain_area > baseDieA.remain_area*3/10) {
            result.indieA_map[cell.name] = true;
            result.dieA.remain_area -= cell.areaA;
        }
        else {
            result.indieA_map[cell.name] = false;
            result.dieB.remain_area -= cell.areaB;
        }
    }
    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result = correct_the_partition(result, cellInfos, libCells, nets);

    result.cost = (result.dieA.remain_area > 0 && result.dieB.remain_area > 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_37(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;
    for (const auto& cell : cellInfos) {
        if (result.dieB.remain_area > baseDieB.remain_area*3/10) {
            result.indieA_map[cell.name] = false;
            result.dieB.remain_area -= cell.areaB;
        }
        else {
            result.indieA_map[cell.name] = true;
            result.dieA.remain_area -= cell.areaA;
        }
    }
    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result = correct_the_partition(result, cellInfos, libCells, nets);

    result.cost = (result.dieA.remain_area > 0 && result.dieB.remain_area > 0) ? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}

PartitionResult initialPartition_net(const vector<CellInfo>& cellInfos,
    const Die baseDieA, const Die baseDieB,
    const pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    const unordered_map<string, vector<connectedCell>>& nets) {

    PartitionResult result;
    result.dieA = baseDieA;
    result.dieB = baseDieB;

    vector<pair<CellInfo, double>> cellWithWeight;
    for (const auto& cell : cellInfos) {
        double totalWeight = 0;
        if (nets.find(cell.name) != nets.end()) {
            for (const auto& conn : nets.at(cell.name)) {
                totalWeight += conn.weight;
            }
        }
        cellWithWeight.push_back({ cell, totalWeight });
    }
    sort(cellWithWeight.begin(), cellWithWeight.end(), [](const pair<CellInfo, double>& a,
        const pair<CellInfo, double>& b) {
            return a.second > b.second;
        });

    for (const auto& cellPair : cellWithWeight) {
        const CellInfo& cell = cellPair.first;
        double scoreA = 0, scoreB = 0;
        if (nets.find(cell.name) != nets.end()) {
            for (const auto& conn : nets.at(cell.name)) {
                if (result.indieA_map.find(conn.name) != result.indieA_map.end()) {
                    if (result.indieA_map.at(conn.name))
                        scoreA += conn.weight;
                    else
                        scoreB += conn.weight;
                }
            }
        }
        bool assignToA = false;
        bool canAssignA = (result.dieA.remain_area >= cell.areaA);
        bool canAssignB = (result.dieB.remain_area >= cell.areaB);
        if (scoreA > scoreB) {
            if (canAssignA)
                assignToA = true;
            else if (canAssignB)
                assignToA = false;
            else
                assignToA = true; // forced assignment
        }
        else if (scoreB > scoreA) {
            if (canAssignB)
                assignToA = false;
            else if (canAssignA)
                assignToA = true;
            else
                assignToA = false;
        }
        else { 
            if (canAssignA && canAssignB) {
                assignToA = (result.dieA.remain_area >= result.dieB.remain_area);
            }
            else if (canAssignA)
                assignToA = true;
            else if (canAssignB)
                assignToA = false;
            else
                assignToA = true;
        }

        if (assignToA) {
            result.indieA_map[cell.name] = true;
            result.dieA.remain_area -= cell.areaA;
        }
        else {
            result.indieA_map[cell.name] = false;
            result.dieB.remain_area -= cell.areaB;
        }
    }

    for (const auto& cell : cellInfos) {
        if (result.indieA_map[cell.name])
            result.cellsA.push_back(cell.name);
        else
            result.cellsB.push_back(cell.name);
    }

    result = correct_the_partition(result, cellInfos, libCells, nets);
    result.cost = (result.dieA.remain_area >= 0 && result.dieB.remain_area >= 0)? computeInitialCost(result.indieA_map, nets) : MAX_double;
    return result;
}



double simulatedAnnealing(double initial_cost, unordered_map<string, vector<connectedCell>>& nets,
    unordered_map<string, bool>& indieA_map, Die& dieA, Die& dieB,
    pair<unordered_map<string, long long>, unordered_map<string, long long>>& libCells,
    vector<string>& cellsA, vector<string>& cellsB) {

    double temperature = 200.0;
    double coolingRate = 0.9999;
    double minTemperature = 0.0001;
    double cost = initial_cost;

    while (temperature > minTemperature) {
        if (cellsA.empty() || cellsB.empty()) break;

        long long selected_idxA, selected_idxB;
        string cellA, cellB;
        bool validSwap = false;
        int failedcnt = 0;
        while (!validSwap && failedcnt < 100000) {
            failedcnt++;
            selected_idxA = rand() % cellsA.size();
            selected_idxB = rand() % cellsB.size();
            cellA = cellsA[selected_idxA];
            cellB = cellsB[selected_idxB];

            if (dieA.remain_area - libCells.first[cellB] + libCells.first[cellA] >= 0 &&
                dieB.remain_area + libCells.second[cellB] - libCells.second[cellA] >= 0) {
                validSwap = true;
            }
        }

        if (failedcnt == 100000) {
            break;
        }
        else
            failedcnt = 0;

        double newCost = cost;

        for (const auto& conn : nets[cellA]) {
            if (conn.name != cellB) {
                if (indieA_map[conn.name] != indieA_map[cellA]) {
                    newCost -= conn.weight;
                }
                else
                    newCost += conn.weight;
            }
        }

        for (const auto& conn : nets[cellB]) {
            if (conn.name != cellA) {
                if (indieA_map[conn.name] != indieA_map[cellB]) {
                    newCost -= conn.weight;
                }
                else
                    newCost += conn.weight;
            }
        }

        double acceptanceProbability = exp((cost - newCost) / temperature);

        if (newCost < cost || ((double)rand() / RAND_MAX) < acceptanceProbability) {
            cost = newCost;
            dieA.remain_area += libCells.first[cellA] - libCells.first[cellB];
            dieB.remain_area += libCells.second[cellB] - libCells.second[cellA];
            cellsA[selected_idxA] = cellB;
            cellsB[selected_idxB] = cellA;
            indieA_map[cellA] = !indieA_map[cellA];
            indieA_map[cellB] = !indieA_map[cellB];
        }

        temperature *= coolingRate;
    }

    return cost;
}


void outputfile(string outfilename, unordered_map<string, vector<connectedCell>> nets, vector<string>cellsA,
    vector<string>cellsB, unordered_map<string, bool>indieA_map, realNet realnet) {

    ofstream outfile(outfilename);

    long long cutsize = 0;
    for (string cellA : cellsA) {
        for (connectedCell conn : nets[cellA]) {
            if (!indieA_map[conn.name]) {
                for (long long i = 0; i < realnet.connected_by_Nets[cellA][conn.name].size(); i++) {
                    if (!realnet.Net_beenUsed[realnet.connected_by_Nets[cellA][conn.name][i]]) {
                        cutsize += realnet.Net_realWeight[realnet.connected_by_Nets[cellA][conn.name][i]];
                        realnet.Net_beenUsed[realnet.connected_by_Nets[cellA][conn.name][i]] = true;
                    }
                }
            }
        }
    }

    outfile << "CutSize " << cutsize << endl;
    outfile << "DieA " << cellsA.size() << endl;
    for (const string& cell : cellsA) {
        outfile << cell << endl;
    }
    outfile << "DieB " << cellsB.size() << endl;
    for (const string& cell : cellsB) {
        outfile << cell << endl;
    }
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: ./hw2 <input.txt> <output.out>" << endl;
        return 1;
    }

    string filename = argv[1];
    string outputfilename = argv[2];
    srand(6);
    gStartTime = chrono::steady_clock::now();


    pair<unordered_map<string, long long>, unordered_map<string, long long>> libCells;
    Die baseDieA, baseDieB;
    unordered_map<string, vector<connectedCell>> nets;
    realNet realnet;
    vector<CellInfo> cellInfos;

    readfile(filename, libCells, baseDieA, baseDieB, nets, realnet, cellInfos);

    PartitionResult partNet = initialPartition_net(cellInfos, baseDieA, baseDieB, libCells, nets);
    gPartNet = partNet;

    PartitionResult partMustok = initialPartition_mustok(cellInfos, baseDieA, baseDieB, libCells, nets);
    //cout << "Mustok strategy initial cost: " << partMustok.cost << endl;
    PartitionResult partSortA   = initialPartition_sortA(cellInfos, baseDieA, baseDieB, libCells, nets);
    //cout << "SortA strategy initial cost: " << partSortA.cost << endl;
    PartitionResult partRandom  = initialPartition_random(cellInfos, baseDieA, baseDieB, libCells, nets);
    //cout << "Random strategy initial cost: " << partRandom.cost << endl;
    PartitionResult partGreedy = initialPartition_greedy(cellInfos, baseDieA, baseDieB, libCells, nets);
    //cout << "Greedy strategy initial cost: " << partGreedy.cost << endl;
    PartitionResult partOriginal = initialPartition_original(cellInfos, baseDieA, baseDieB, libCells, nets);
    //cout << "Original strategy initial cost: " << partOriginal.cost << endl;
	PartitionResult part73 = initialPartition_73(cellInfos, baseDieA, baseDieB, libCells, nets);
	PartitionResult part37 = initialPartition_37(cellInfos, baseDieA, baseDieB, libCells, nets);

    if(timeoutReached()){
        outputfile(outputfilename, nets, partNet.cellsA, partNet.cellsB, partNet.indieA_map, realnet);
        return 0;
    }


    vector<PartitionResult> partitions = { partNet, partMustok, partSortA, partRandom, partGreedy, partOriginal, part73, part37 };

    double bestFinalCost = MAX_double;
    vector<string> bestCellsA;
    vector<string> bestCellsB;
    unordered_map<string, bool> bestIndieA_map;
    Die bestDieA, bestDieB;

	for(int i = 0; i<3; i++){
		for (const auto& part : partitions) {
			if (part.cost == MAX_double)
				continue; 

			vector<string> tmpCellsA = part.cellsA;
			vector<string> tmpCellsB = part.cellsB;
			unordered_map<string, bool> tmpIndieA_map = part.indieA_map;
			Die tmpDieA = part.dieA;
			Die tmpDieB = part.dieB;
			double simCost = simulatedAnnealing(part.cost, nets, tmpIndieA_map, tmpDieA, tmpDieB, libCells, tmpCellsA, tmpCellsB);
			//cout << "Simulated annealing cost: " << simCost << endl;
			if (simCost < bestFinalCost) {
				bestFinalCost = simCost;
				bestCellsA = tmpCellsA;
				bestCellsB = tmpCellsB;
				bestIndieA_map = tmpIndieA_map;
				bestDieA = tmpDieA;
				bestDieB = tmpDieB;
			}
            if(timeoutReached())
                break;
		}
        if(timeoutReached())
            break;
	}

    //cout << "Selected final cost after simulated annealing: " << bestFinalCost << endl;
    
    outputfile(outputfilename, nets, bestCellsA, bestCellsB, bestIndieA_map, realnet);

    return 0;
}
