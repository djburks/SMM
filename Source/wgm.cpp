# include <fstream>
# include <string>
# include <sstream>
# include <iostream>
# include <unordered_map>
# include <vector>
# include <math.h>
# include <algorithm>

//int t = 0;

std::unordered_map<char,int> nucloc = {
    {'A',1},
    {'T',2},
    {'G',3},
    {'C',4}
};

std::unordered_map<char,int> rnucloc = {
    {'A',2},
    {'T',1},
    {'G',4},
    {'C',3}
};

struct kmer {
    double counts[5] = {4,1,1,1,1};
    double * total = &counts[0];
};
    
struct kmer kset[16777216];

int index(std::string k) {
    int idex = 0;
    int len = k.length();
    for(int i=0; i<=len; i++) {
        idex = idex + ((pow(4,len-i-1))*nucloc[k[i]]);
    }
    return idex;
}

int rindex(std::string k) {
    int ridex = 0;
    int len = k.length();
    for(int i=0; i<=len; i++) {
        ridex = ridex + (pow(4,len-i-1)*rnucloc[k[i]]);
    }
    return ridex;
}

int main(int argc, char *argv[]) {
    

// Crude Fasta Import
    int order = std::stoi(argv[2]);
    
    std::string line;
    std::string genome = "";
    std::ifstream infile(argv[1]);
    while (std::getline(infile,line)) {
        if (line[0] != '>') {
            genome.append(line);
        }
    }
    
    

// Iterate through genome, build counts, and convert to probabilities (the model.).
    

    int giter = genome.length() - order;
    std::string r(order,'A');
    std::string e(order,'C');
    
    int reducer = index(r);
    int ender = index(e);
    long int totalc = giter + (pow(4,order)*4);

    for(int i=0; i<giter; i++) {
        std::string kmer =  genome.substr(i,order);
        char tnuc = genome[i+order];
        int idex = (index(kmer) - reducer);
        kset[idex].counts[nucloc[tnuc]] = kset[idex].counts[nucloc[tnuc]] + 1;
        kset[idex].counts[0]++;

    }
    
    

    for(int i=0; i<=(ender-reducer); i++) {
        kset[i].counts[4] = log10(kset[i].counts[4]/kset[i].counts[0]);
        kset[i].counts[3] = log10(kset[i].counts[3]/kset[i].counts[0]);
        kset[i].counts[2] = log10(kset[i].counts[2]/kset[i].counts[0]);
        kset[i].counts[1] = log10(kset[i].counts[1]/kset[i].counts[0]);
        kset[i].counts[0] = log10(kset[i].counts[0]/totalc);
    }
    


    
        
// Take in fasta read file, store in ordered map.

    std::ifstream readfile(argv[3]);
    std::string readline;
    std::string curentry,probline;
    probline.append(argv[1]);
    std::vector<std::string> metafasta;
    
    while (std::getline(readfile,readline)) {
        if (readline[0] == '>') {
            metafasta.push_back("");
        }
        else {
            (metafasta[metafasta.size() - 1]).append(readline);
        }
    }
    
// Cycle through metafasta map and score reads.  Take best value (reverse or forward).
    
    double forprob,revprob;
    
    for(auto& it: metafasta) {
        std::string fread = it;
        std::string rread = it;
        std::reverse(rread.begin(),rread.end());
        std::string&& fchunk = fread.substr(0,order);
        std::string&& rchunk = rread.substr(0,order);
        
        int fidex = 0;
        int ridex = 0;
        for(int j=0; j<order; j++) {
                fidex = fidex + (pow(4,order-j-1)*nucloc[fchunk[j]]);
                ridex = ridex + (pow(4,order-j-1)*rnucloc[rchunk[j]]);
        }
        double forprob = kset[fidex-reducer].counts[0];
        double revprob = kset[ridex-reducer].counts[0];
        
// Cycle through read and compute running score
        for(int i=0;i<fread.length()-order;i++) {
            std::string&& fchunk = fread.substr(i,order);
            std::string&& rchunk = rread.substr(i,order);

            
            int fidex = 0;
            int ridex = 0;
            for(int j=0; j<order; j++) {
                fidex = fidex + (pow(4,order-j-1)*nucloc[fchunk[j]]);
                ridex = ridex + (pow(4,order-j-1)*rnucloc[rchunk[j]]);
            }
            fidex = fidex - reducer;
            ridex = ridex - reducer;
            //forprob = forprob + kset[index(fread.substr(i,order))-reducer].counts[nucloc[fread[order]]];
            forprob = forprob + kset[fidex].counts[nucloc[fread[i+order]]];
            revprob = revprob + kset[ridex].counts[rnucloc[rread[i+order]]];
        }
        
        probline.append("\t" + std::to_string(std::max(forprob,revprob)));
       // t++;
       // std::cout << t << std::endl;

        
        
    }
    
    std::cout << probline << std::endl;

    return 0;
    
}
