#include "../googletest-main/googletest/include/gtest/gtest.h"
#include "../include/cluster.h"
#include <stdio.h>
#include <math.h>

namespace {
    class ClusterTest : public ::testing::Test {
        protected:

        ClusterTest() {
            vector<double> memorySizes = Cluster::buildNLevelMemorySizes(mems, processorCounts);
            clusterFixture ->setMemorySizes(memorySizes);
            Cluster::setFixedCluster(clusterFixture);
        }

        ~ClusterTest() override {
            //delete clusterFixture;

        }
    unsigned int num_processors = 1500;
    Cluster *clusterFixture = new Cluster(num_processors,false);
    vector<unsigned int> processorCounts = {100,500,900};
    vector<double> mems = {200,100,10};
    };

    bool double_equal(double d1, double d2){
        if (abs(d1-d2)<1.0e-10) return true;
        return false;
    }

}


TEST_F(ClusterTest, buildStatic2LevelClusterIsHeterogeneous){
    Cluster::buildStatic2LevelCluster(100,10);
    Cluster *cluster = Cluster::getFixedCluster();
    cluster->printInfo();
    EXPECT_FALSE(cluster->isHomogeneous());
    
}

TEST_F(ClusterTest, buildStatic2LevelClusterCorrectSize){
    Cluster::buildStatic2LevelCluster(100,10);
    Cluster *cluster = Cluster::getFixedCluster();

    EXPECT_EQ(cluster->getNumberProcessors(), 3758);
}

TEST_F(ClusterTest, buildStatic2LevelClusterIs2Level){
    Cluster::buildStatic2LevelCluster(100,10);
    Cluster *cluster = Cluster::getFixedCluster();
    vector<Processor *> proc = cluster->getProcessors();

    double mem1 = proc.front()->getMemorySize();
    double mem2 = proc.back()->getMemorySize();
    for (Processor *p:proc){
        double mem = p->getMemorySize();
        EXPECT_TRUE(double_equal(mem,mem1)||double_equal(mem,mem2));
    }
}

TEST_F(ClusterTest, buildStatic3LevelClusterIsHeterogeneous){
    Cluster::buildStatic3LevelCluster(100,10);
    Cluster *cluster = Cluster::getFixedCluster();
    cluster->printInfo();
    EXPECT_FALSE(cluster->isHomogeneous());
    
}

TEST_F(ClusterTest, buildStatic3LevelClusterCorrectSize){
    Cluster::buildStatic3LevelCluster(100,10);
    Cluster *cluster = Cluster::getFixedCluster();

    EXPECT_EQ(cluster->getNumberProcessors(), 1475);
}

TEST_F(ClusterTest, buildStatic3LevelClusterIs3Level){
    Cluster::buildStatic2LevelCluster(100,10);
    Cluster *cluster = Cluster::getFixedCluster();
    vector<Processor *> proc = cluster->getProcessors();

    double mem1 = proc.front()->getMemorySize();
    double mem2 = proc.back()->getMemorySize();
    double mem3 = 0;
    for (Processor *p:proc){
        double mem = p->getMemorySize();
        if ( !double_equal(mem1,mem) && !double_equal(mem2,mem)){
            mem3 = mem;
            break;
        }
    }
    for (Processor *p:proc){
        double mem = p->getMemorySize();
        EXPECT_TRUE(double_equal(mem,mem1)||double_equal(mem,mem2)||double_equal(mem,mem3));
    }
}

TEST_F(ClusterTest, testBuildNLevelMemorySizescorrectSize){
    Cluster * cluster = Cluster::getFixedCluster();
    EXPECT_EQ(num_processors, cluster->getNumberProcessors());
}

TEST_F(ClusterTest, testBuildNLevelIs3Level){
    Cluster * cluster = Cluster::getFixedCluster();
    vector<Processor *> proc = cluster->getProcessors();

    double mem1 = proc.front()->getMemorySize();
    double mem2 = proc.back()->getMemorySize();
    double mem3 = 0;
    for (Processor *p:proc){
        double mem = p->getMemorySize();
        if ( !double_equal(mem1,mem) && !double_equal(mem2,mem)){
            mem3 = mem;
            break;
        }
    }
    for (Processor *p:proc){
        double mem = p->getMemorySize();
        EXPECT_TRUE(double_equal(mem,mem1)||double_equal(mem,mem2)||double_equal(mem,mem3));
    }
}

TEST_F(ClusterTest, testBuildNLevelMemorySizesCorrectlyAllocated){
    Cluster * cluster = Cluster::getFixedCluster();
    vector<Processor *> proc = cluster->getProcessors();
    const int size = mems.size();
    vector<int> cnts(size, 0);
    for (Processor *p:proc){
        double memory = p->getMemorySize();
        for (int i  =0; i<size; i++){
            if (double_equal(mems.at(i), memory)){
                cnts[i]++;
                break;
            }
        } 
    }
    for (int i =0; i<size; i++){
        EXPECT_EQ(cnts[i],processorCounts[i]);
    }
    
}

TEST_F(ClusterTest, estBuildNLevelMemorySizesIsHeterogeneous){
    Cluster *cluster = Cluster::getFixedCluster();
    cluster->printInfo();
    EXPECT_FALSE(cluster->isHomogeneous());
    
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
