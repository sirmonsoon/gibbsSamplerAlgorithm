#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <climits>
#include <numeric>


using namespace std;

//Given Motifs we can construct Profile(Motifs)
//Using pseudocounts to avoid 0 probability option
//Exclude the i-th motif
vector<vector<double>> createProfile(const vector<string>& motifs, int k, int exclude = -1){
    vector<vector<double>> profile(4, vector<double>(k, 1.0));
    int t = motifs.size();

    for(int i = 0; i < t; i++){
        if (i == exclude) continue;
        for(int j = 0; j < k; j++){
            char nucleotide = motifs[i][j];
            if(nucleotide == 'A') 
            profile[0][j]++;
            else if(nucleotide == 'C')
            profile[1][j]++;
            else if(nucleotide == 'G')
            profile[2][j]++;
            else if(nucleotide == 'T')
            profile[3][j]++;
        }
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < k; j++){
            profile[i][j] /= (t - (exclude >= 0 ? 1 : 0) + 4.0);
        }
    }
    return profile;
}

//Scoring Motifs based on consensus.
int scoreMotifs(const vector<string>& motifs, int k){
    int score = 0;
    for(int j = 0; j < k; j++){
        vector<int> count(4, 0);
        for(int i = 0; i < motifs.size(); i++){
            char nucleotide = motifs[i][j];
            if(nucleotide == 'A')
            count[0]++;
            else if(nucleotide == 'C')
            count[1]++;
            else if(nucleotide == 'G')
            count[2]++;
            else if(nucleotide == 'T')
            count[3]++;
        }
        score += motifs.size() - *max_element(count.begin(), count.end());
    }
    return score;
}

//Randomly select kmers from the given sequence 
vector<string> randomlySelectedKmers(const vector<string>& dna, int k){
    vector<string> motifs;
    for(const string& seq : dna){
        int randIndex = rand() % (seq.size() - k + 1);
        motifs.push_back(seq.substr(randIndex, k));
    }
    return motifs;
}

//Sample kmer for a sequence using given profile
string sampleKmer(const string& sequence, const vector<vector<double>>& profile, int k){
    vector<double> probs(sequence.size() - k + 1, 0.0);

    for(int i = 0; i <= sequence.size() - k; i++){
        string kmer = sequence.substr(i, k);
        double prob = 1.0;
        for(int j = 0; j < k; j++){
            char nucleotide = kmer[j];
            if(nucleotide == 'A')
                prob *= profile[0][j];
            else if(nucleotide == 'C')
                prob *= profile[1][j];
            else if(nucleotide == 'G')
                prob *= profile[2][j];
            else if(nucleotide == 'T')
                prob *= profile[3][j];
            }
            probs[i] = prob;
    }

    double totalProb = accumulate(probs.begin(), probs.end(), 0.0);
    for(double& p : probs){
        p /= totalProb;
    }
    double randProb = ((double) rand() / RAND_MAX);
    double cumulativeProb = 0.0;
    for(int i = 0; i < probs.size(); i++){
        cumulativeProb += probs[i];
        if(randProb <= cumulativeProb){
            return sequence.substr(i, k);
        }
    }
    return sequence.substr(0, k);
}

//Iteratively improves motifs based on profile matrix
vector<string> gibbsSampler(const vector<string>& dna, int k, int t, int r){
    vector<string> motifs = randomlySelectedKmers(dna, k);
    vector<string> bestMotifs = motifs;
    int bestScore = scoreMotifs(motifs, k);

    //For j <- 1 to N
    for(int iterations = 0; iterations < r; iterations++){
        //i <- RANDOM(t)
        int i = rand() % t;
        //Profile <- profile matrix formed from all strings in Motifs except for Motif[i]
        vector<vector<double>> profile = createProfile(motifs, k, i);
        //Motif[i] <- Profile-randomly generated k-mer in the i-th sequence.
        motifs[i] = sampleKmer(dna[i], profile, k);

        int currScore = scoreMotifs(motifs, k);
        //If score(motifs) < score(bestmotifs) then bestmotifs <- motifs
        if(currScore < bestScore){
            bestMotifs = motifs;
            bestScore = currScore;
        }
    }
    return bestMotifs;
}

int main(int argc, char* argv[])
{
    srand(time(0)); 
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << "<input_file>" << endl;
        return 1;
    }

    string inputFileName = argv[1];
    ifstream inputFile(inputFileName);
    if (!inputFile.is_open()){
        cerr << "File Read Error" << endl;
        return 1;
    }

    int k, t, r;
    inputFile >> k >> t >> r;
    vector<string> dna(t);
    for(int i = 0; i < t; i++){
        inputFile >> dna[i];
    }
    inputFile.close();

    vector<string> bestMotifs; //gibbsSampler(dna, k, t, r);
    int bestScore = INT_MAX;

    //Running algorithm
    for(int i = 0; i < 30; i++){
        vector<string> motifs = gibbsSampler(dna, k, t, r);
        int currScore = scoreMotifs(motifs, k);

        if(currScore < bestScore){
            bestMotifs = motifs;
            bestScore = currScore;
        }
    }
    
    //Output file name extraction from input file.
    string qNum = "2";
    string testCaseNum = inputFileName.substr(inputFileName.find('_') + 1, inputFileName.find('.') - inputFileName.find('_') - 1);
    stringstream outputFileName;

    outputFileName << "sol_q" << qNum << "_t" << testCaseNum;

    //Write output to output file.
    ofstream outputFile(outputFileName.str());
    if(!outputFile){
        cerr << "Error: Unable to open output file " << outputFileName.str() << endl;
        return 1;
    }

    for(const string& motif : bestMotifs){
        outputFile << motif << endl;
        //Testing output
        //cout << motif << endl;
    }
    outputFile.close();

    return 0;

}