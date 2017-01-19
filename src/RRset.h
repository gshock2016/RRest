#include <fstream>

class Math{
    public:
        static double log2(int n){
            return log(n) / log(2);
        }
        static double logcnk(int n, int k) {
            double ans = 0;
            for (int i = n - k + 1; i <= n; i++)
            {
                ans += log(i);
            }
            for (int i = 1; i <= k; i++)
            {
                ans -= log(i);
            }
            return ans;
        }
};
class RRset
{
    private:
        //static InfGraph g;
        //static int k;
        //static map<string, string> arg;

        static double step1(InfGraph &g, const Argument & arg)
        {
            double epsilon_prime = arg.epsilon * sqrt(2);
            Timer t(1, "step1");
            for (int x = 1; ; x++)
            {
                int64 ci = (2+2/3 * epsilon_prime)* ( log(g.n) + Math::logcnk(g.n, arg.k) + log(Math::log2(g.n))) * pow(2.0, x) / (epsilon_prime* epsilon_prime);
                INFO(ci);
                g.build_hyper_graph_r(ci, arg);
                g.build_seedset(arg.k);
                double ept = g.InfluenceHyperGraph()/g.n;
                double estimate_influence = ept * g.n;
                INFO(x, estimate_influence);
                cout<<"ept = "<<ept<<endl<<"1/2^x = "<<1 / pow(2.0, x)<<endl;
                if (ept > 1 / pow(2.0, x))
                {
                    double OPT_prime = ept * g.n / (1+epsilon_prime);
                    //INFO("step1", OPT_prime);
                    //INFO("step1", OPT_prime * (1+epsilon_prime));
                    return OPT_prime;

                }
            }
            ASSERT(false);
            return -1;
        }
        static double step2(InfGraph &g, const Argument & arg, double OPT_prime)
        {
            Timer t(2, "step2");
            ASSERT(OPT_prime > 0);
            double e = exp(1);
            double alpha = sqrt(log(g.n) + log(2));
            double beta = sqrt((1-1/e) * (Math::logcnk(g.n, arg.k) + log(g.n) + log(2)));

            int64 R = 2.0 * g.n *  sqr((1-1/e) * alpha + beta) /  OPT_prime / arg.epsilon / arg.epsilon ;
            cout<<"R = "<<R<<endl;
            //R/=100;
            g.build_hyper_graph_r(R, arg);
            g.build_seedset(arg.k);
            double opt = g.InfluenceHyperGraph();
            return opt;
        }
    public:
        static void InfluenceMaximize(InfGraph &g, const Argument &arg)
        {
            Timer t(100, "InfluenceMaximize(Total Time)");

            INFO("########## Step1 ##########");

            // debugging mode lalala
            double OPT_prime;
            OPT_prime = step1(g, arg ); //estimate OPT



            INFO("########## Step2 ##########");


            double opt_lower_bound = OPT_prime;
            INFO(opt_lower_bound);
            step2(g, arg, OPT_prime);
            INFO("step2 finish");

        }
    /*
    * Output RRsets
    */
        static void getRRsets(InfGraph &g, const Argument &arg, int threshold){

            double OPT_prime = step1(g, arg);

            double e = exp(1);
            double alpha = sqrt(log(g.n) + log(2));
            double beta = sqrt((1-1/e) * (Math::logcnk(g.n, arg.k) + log(g.n) + log(2)));

            int64 R = 2.0 * g.n *  sqr((1-1/e) * alpha + beta) /  OPT_prime / arg.epsilon / arg.epsilon ;

            cout<<"R= "<<R<<endl;

            g.build_hyper_graph_r(R, arg);
            ofstream myfile;
            myfile.open ("./output/rrset.csv");

            int nodeID;

            vector<int> infValues;
            for(int i =0; i < g.n; i++){
                infValues.push_back(0);
            }
/*
 * RRset ID start from 0
 * Output inlfuence matrix, infadjlist, rrsets
 * Example:
 * RRset:
 * -------------------------------------------------------
 * seed node | sequence of influenced node (BFS sequence)
 *     2     | 1, 3,
 *     1     | 2, 3
 *     3     | 2,
 * -------------------------------------------------------
 * InfValue:
 * -------------------------------------------------------
 * node ID   | influence value (number of rrsets it belongs to)
 *     1     | 2
 *     2     | 3
 *     3     | 3
 * -------------------------------------------------------
 * InfAdjList:
 * -------------------------------------------------------
 * node ID   | rrsets the node belongs to
 *     1     | 1, 2
 *     2     | 1, 2, 3
 *     3     | 1, 2, 3
 * -------------------------------------------------------
 * InfMatrix
 * -------------------------------------------------------
 * node ID   | if node i belong to rrset j, then (i, j) = 1, otherwise (i,j)=0
 *     1     | 1, 1, 0
 *     2     | 1, 1, 1
 *     3     | 1, 1, 1
 * -------------------------------------------------------
 */
            for (int i = 0; i < (int) g.hyperGT.size();i++){
                int rrsetID = i;
                for(int j = 0; j < (int) g.hyperGT[i].size(); j++){
//                    cout<<g.hyperGT[i][j]<<", ";
                    nodeID = g.hyperGT[i][j];
                    myfile <<nodeID<<",";
                    g.infAdjList[nodeID].push_back(rrsetID);
                    g.infmatrix[nodeID][rrsetID]=1;
                    infValues[nodeID] = infValues[nodeID] + 1;
                }
                myfile << "\n";
//                cout<<endl;
            }
            myfile.close();

            ofstream myfile2;
            myfile2.open("./output/infvalue.csv");
            ofstream myfile3;
            myfile3.open("./output/infmatrix.csv");
            ofstream myfile4;
            myfile4.open("./output/infadjlist.csv");
            for (int i = 0; i < g.n; i++){
//                myfile2<<i<<","<<infValues[i]<< "\n";
                if(infValues[i] > threshold) {
                    myfile2<<i<<","<<infValues[i]<< "\n";
                    for (int j = 0; j < (int) g.infmatrix[i].size(); j++) {
                        myfile3 << g.infmatrix[i][j] << ",";
                    }
                    for (int j = 0; j < g.infAdjList[i].size(); j++) {
//                        cout<<(int) g.infAdjList[i].size()<<endl;
                        myfile4 << g.infAdjList[i][j] << ",";
                    }
                    myfile3 << "\n";
                    myfile4 << "\n";
                }

            }
            myfile2.close();
            myfile3.close();
            myfile4.close();
        }

    static void getRRsets(InfGraph &g, int R, const Argument &arg, int threshold){

        g.build_hyper_graph_r(R, arg);
        ofstream myfile;
        myfile.open ("./output/rrset.csv");

        int nodeID;

        vector<int> infValues;
        for(int i =0; i < g.n; i++){
            infValues.push_back(0);
        }

        for (int i = 0; i < (int) g.hyperGT.size();i++){
            int rrsetID = i;
            for(int j = 0; j < (int) g.hyperGT[i].size(); j++){
//                    cout<<g.hyperGT[i][j]<<", ";
                nodeID = g.hyperGT[i][j];
                myfile <<nodeID<<",";
                g.infAdjList[nodeID].push_back(rrsetID);
                g.infmatrix[nodeID][rrsetID]=1;
                infValues[nodeID] = infValues[nodeID] + 1;
            }
            myfile << "\n";
//                cout<<endl;
        }
        myfile.close();

        ofstream myfile2;
        myfile2.open("./output/infvalue.csv");
        ofstream myfile3;
        myfile3.open("./output/infmatrix.csv");
        ofstream myfile4;
        myfile4.open("./output/infadjlist.csv");
        for (int i = 0; i < g.n; i++){
//                myfile2<<i<<","<<infValues[i]<< "\n";
            if(infValues[i] > threshold) {
                myfile2<<i<<","<<infValues[i]<< "\n";
//                for (int j = 0; j < (int) g.infmatrix[i].size(); j++) {
//                    myfile3 << g.infmatrix[i][j] << ",";
//                }
                for (int j = 0; j < g.infAdjList[i].size(); j++) {
//                        cout<<(int) g.infAdjList[i].size()<<endl;
                    myfile4 << g.infAdjList[i][j] << ",";
                }
                myfile3 << "\n";
                myfile4 << "\n";
            }

        }
        myfile2.close();
        myfile3.close();
        myfile4.close();
    }

};

