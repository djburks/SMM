# include <fstream>
# include <string>
# include <sstream>
# include <iostream>
# include <unordered_map>
# include <vector>
# include <math.h>
# include <algorithm>


// Transition nucleotide index.

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

// Fasta import and filter function.
// Skips id lines, returns only standard ATGC nucleotides in uppercase.

std::string fasta(std::string infile) {
    std::ifstream fna(infile);
    std::string line;
    std::string genome;
    while(std::getline(fna,line)) {
        if (line[0] == '>') {
            continue;
        }
        for(char const &c: line) {
            if (c == 'A' || c == 'a') {genome+=('A');}
            if (c == 'T' || c == 't') {genome+=('T');}
            if (c == 'G' || c == 'g') {genome+=('G');}
            if (c == 'C' || c == 'c') {genome+=('C');}
        }
        
        
    }
    return genome;
}

// Kmer index functions.

int index(std::string k) {
    int idex = 0;
    int len = k.length();
    for(int i=0; i<=len; i++) {
        idex = idex + ((pow(4,len-i-1))*nucloc[k[i]]);
    }
    return idex;
}

// Build multidimensional array model for storage of counts.
// 1-4 for transitional specific counts, 0 for total.

std::vector<std::vector<double>> smm(std::string genomefile, int order) {
    int klim = (int)pow(4,order) + 1;
    std::vector<std::vector<double>> model(5,std::vector<double> (klim,1));
    std::string line;
    std::string genome = "";
    std::ifstream infile(genomefile);
    while (std::getline(infile,line)) {
        if (line[0] == '>') {
            continue;
        }
        for(char const &c: line) {
            if (c == 'A' || c == 'a') {genome+=('A');}
            if (c == 'T' || c == 't') {genome+=('T');}
            if (c == 'G' || c == 'g') {genome+=('G');}
            if (c == 'C' || c == 'c') {genome+=('C');}
        }
    }

    // Iterate through genome, build counts, and generate probabilities.
    int giter = genome.length() - order;
    std::string r(order,'A');
    std::string e(order,'C');
    int reducer = index(r);
    int ender = index(e);
    int totalc = giter + (pow(4,order)*4);

    
    for(int i=0; i<giter; i++) {
        std::string&& kmer =  genome.substr(i,order);
        char tnuc = genome[i+order];
        int idex = (index(kmer) - reducer);
        model[nucloc[tnuc]][idex]++;
        model[0][idex]++;
    }
    

    for(int i=0; i<=(ender-reducer); i++) {
        model[4][i] = log10((model[4][i])/(model[0][i] + 3));
        model[3][i] = log10((model[3][i])/(model[0][i] + 3));
        model[2][i] = log10((model[2][i])/(model[0][i] + 3));
        model[1][i] = log10((model[1][i])/(model[0][i] + 3));
        model[0][i] = log10((model[0][i] + 3)/totalc);
    }

    return model;
}

// Read in metafasta file.
// Store filtered reads in vector.

std::vector<std::string> readlibrary(std::string infile) {
    std::vector<std::string> metafasta;
    std::ifstream readfile(infile);
    std::string readline,filterline;
    while (std::getline(readfile,readline)) {
        if (readline[0] == '>') {
            metafasta.push_back("");
            filterline = "";
        }
        else {
            filterline = "";
            for(char const &c: readline) {
                if (c == 'A' || c == 'a') {filterline+=('A');}
                if (c == 'T' || c == 't') {filterline+=('T');}
                if (c == 'G' || c == 'g') {filterline+=('G');}
                if (c == 'C' || c == 'c') {filterline+=('C');}
            }
            (metafasta[metafasta.size() - 1]).append(filterline);
        }
    }
    return metafasta;
}



int rindex(std::string k) {
    int ridex = 0;
    int len = k.length();
    for(int i=0; i<=len; i++) {
        ridex = ridex + (pow(4,len-i-1)*rnucloc[k[i]]);
    }
    return ridex;
}

// Normalization for average 12th order scores based on read length across all models.

double normScore(float rawscore,int readlen) {
    float denom;
    denom = (-0.605185747917117*readlen) + 0.228066296087579;
    return rawscore/denom;
}

int main(int argc, char *argv[]) {
    std::string readfile = argv[2];
    int order = std::stoi(argv[3]);
    std::string scoretype = argv[4];
    std::vector<std::string> mfasta = readlibrary(readfile);
    std::string r(order,'A');
    std::string e(order,'C');
    int reducer = index(r);
    int ender = index(e);
    // Store index vectors for metafasta for rapid-lookup.
    std::vector<std::vector<int>> row;
    std::vector<std::vector<int>> column;
    std::vector<std::vector<int>> rcolumn;
    std::vector<std::vector<int>> rrow;
    std::vector<int> mflen;
    
    for (auto& it: mfasta) {
        std::vector<int> temprow,tempcol,temprrow,temprcol;
        temprow.push_back(0);
        temprrow.push_back(0);
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
        fidex = fidex - reducer;
        ridex = ridex - reducer;
        tempcol.push_back(fidex);
        temprcol.push_back(ridex);
        
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
            
            temprow.push_back(nucloc[fread[i+order]]);
            temprrow.push_back(rnucloc[rread[i+order]]);
            
            tempcol.push_back(fidex);
            temprcol.push_back(ridex);
        
        }
        row.push_back(temprow);
        column.push_back(tempcol);
        rcolumn.push_back(temprcol);
        rrow.push_back(temprrow);
        mflen.push_back(it.length());
    }
    
// For every genome in the linked genome fasta list, cycle through and build model.
// Use fasta lookups to generate score.
    std::vector<std::string> genomes;
    std::string genomeline;
    std::ifstream genomefastafile(argv[1]);
    while(std::getline(genomefastafile,genomeline)) {
        genomes.push_back(genomeline);
    }

    if (scoretype == "norm") {
            for (auto& genome: genomes) {
                std::string probline = genome;
                std::vector<std::vector<double>> smmmodel = smm(genome,order);
                for(int i2=0; i2 < row.size(); i2++) {
                    double forprob = 0;
                    double revprob = 0;
                    for(int i3=0; i3 < row[i2].size(); i3++) {
                        forprob = forprob + smmmodel[row[i2][i3]][column[i2][i3]];
                        revprob = revprob + smmmodel[rrow[i2][i3]][rcolumn[i2][i3]];
                    }
                    probline.append("\t" + std::to_string(normScore(std::max(forprob,revprob),mflen[i2])));
                }
                std::cout << probline << "\n";

            }
    }
    
    else {
            for (auto& genome: genomes) {
                std::string probline = genome;
                std::vector<std::vector<double>> smmmodel = smm(genome,order);
                for(int i2=0; i2 < row.size(); i2++) {
                    double forprob = 0;
                    double revprob = 0;
                    for(int i3=0; i3 < row[i2].size(); i3++) {
                        forprob = forprob + smmmodel[row[i2][i3]][column[i2][i3]];
                        revprob = revprob + smmmodel[rrow[i2][i3]][rcolumn[i2][i3]];
                    }
                    probline.append("\t" + std::to_string(std::max(forprob,revprob)));
                }
                std::cout << probline << "\n";

            }
    }
    

    
    return 0;
}
