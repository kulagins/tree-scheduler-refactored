
#include "../googletest-main/googletest/include/gtest/gtest.h"
#include "../include/cluster.h"
#include "../include/tree.h"
#include <stdio.h>
#include <math.h>
bool double_ge(double v1, double v2){
    return v1-v2 > -1e-10;

}
namespace {
    class DistributeTest : public ::testing::Test {
        protected:

        DistributeTest() {
            vector <unsigned int> processorCounts = {2,2,1};
            vector <double> mems = {10,7,5};
            vector <double> speeds = {5,8,10};

            Cluster *cluster = new Cluster(&processorCounts,&mems,&speeds);
            Cluster::setFixedCluster(cluster);
            

            tree = new Tree();

            Task *t6 = new Task(10, 5, 5, true);
            t6->setId(6);
            Task * t5 = new Task(6,5,5,100.0);
            t5->setId(5);
            Task * t4 = new Task(6,5,5,8.0);
            t4->setId(4);
            Task * t3 = new Task(5,5,5,9.0);
            t3->setId(3);
            Task * t2 = new Task(4,5,5,5.0);
            t2->setId(2);
            Task * t1 = new Task(4,5,5,5.0);
            t1->setId(1);

            tree->addRoot(t6);
            tree->addTask(t5);
            tree->addTask(t4);
            tree->addTask(t3);
            tree->addTask(t2);
            tree->addTask(t1);

            t1->setParent(t4);
            t4->addChild(t1);

            t2->setParent(t4);
            t4->addChild(t2);

            t3->setParent(t5);
            t5->addChild(t3);
            t4->setParent(t6);
            t6->addChild(t4);
            t5->setParent(t6);
            t6->addChild(t5);

            for (auto task: *tree->getTasks()){
                for(auto proc:cluster->getProcessors() ){
                    if(double_ge(proc->getMemorySize(),task->getNodeWeight())){
                        task->addFeasibleProcessor(proc);
                    }
                }
            }

        }

        ~DistributeTest() override {

        }
        Tree *tree;

    };
}

//TEST_F(DistributeTest, )

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
