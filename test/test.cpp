//
// Created by kulaginasv on 22.11.21.
//

#include "../googletest-main/googletest/include/gtest/gtest.h"
#include "../include/tree.h"
#include <stdio.h>
const int DIVISION_LENGTH = 25;
namespace {

    class TreeTest : public ::testing::Test {
        /*
        The tree created here has a long chain that is about 25% of the number of tasks
        another 25% are going to be fork-structured, i.e. there are going to be leaves who's parent is the root
        the last 50% are going to be a well balanced binary tree
        
        */
    protected:

        TreeTest() {
            // initialize all tasks
            Task * task;
            const int NUM_TASKS = 4*DIVISION_LENGTH +1;
            Task * root = new Task(1, 1, 1, true);
            root->setId(1);
            tree->addRoot(root);

            for (int i = 2; i<=NUM_TASKS; i++){
                task = new Task(i, i, i, false);
                task->setId(i);
                tree->addTask(task); 
            }
            
            // build fork part
            // every third edge will be broken
            int curr_task_id = 2;
            for (curr_task_id; curr_task_id<=DIVISION_LENGTH+3; curr_task_id++){
                task = tree->getTask(curr_task_id);
                task->setParent(root);
                root->addChild(task); 

                if (curr_task_id%3 == 0){
                    task->breakEdge();
                }
            }

            // build chain part
            // one edge by arount 80% of the chain will be broken
            Task * parent;
            parent = tree->getTask(2);
            task = tree->getTask(curr_task_id);
            task->setParent(parent);
            parent->addChild(task);  
            curr_task_id ++;

            for (curr_task_id; curr_task_id<(2*DIVISION_LENGTH)+3;curr_task_id++){
                task = tree->getTask(curr_task_id);
                parent =tree->getTask(curr_task_id-1); 
                task->setParent(parent);
                parent->addChild(task); 

                if (curr_task_id == curr_task_id+(int)((curr_task_id<(2*DIVISION_LENGTH)+3 - curr_task_id)*0.8)) {
                    task->breakEdge();
                }
            }

            // build binary-tree part
            // every fourth edge will be broken
            int id_before_binary_tree = curr_task_id -2;
            parent = tree->getTask(DIVISION_LENGTH+3);
            task = tree->getTask(curr_task_id);
            task->setParent(parent);
            parent->addChild(task);  
            curr_task_id ++;

            task = tree->getTask(curr_task_id);
            task->setParent(parent);
            parent->addChild(task);  
            curr_task_id ++;


            for (curr_task_id; curr_task_id<=NUM_TASKS; curr_task_id++){
                parent = tree->getTask((int)((curr_task_id-id_before_binary_tree)/2)+id_before_binary_tree);
                task = tree->getTask(curr_task_id);
                task->setParent(parent);
                parent->addChild(task);  

                if (curr_task_id%4 == 0){
                    task->breakEdge();
                }
            }

        


        }

        ~TreeTest() override {
            delete tree;
        }
    Tree * tree = new Tree();

    };

}

TEST_F(TreeTest, emptyTreeIsEmpty){
    Tree * t = new Tree();
    EXPECT_EQ(t->getSize(),0);
}

TEST_F(TreeTest,testFixtureStructure){
    cout << tree->getSize()<<endl;
    EXPECT_EQ(tree->getSize(),4*DIVISION_LENGTH +1);
    EXPECT_EQ(tree->getRoot()->getChildren()->size(),DIVISION_LENGTH +2);
    cout << tree->getTask((int)1.5*DIVISION_LENGTH)->getId()<<endl;
    EXPECT_EQ(tree->getTask((int)(1.5*DIVISION_LENGTH))->getChildren()->size(),1);
    EXPECT_EQ(tree->getTask((int)(2.5*DIVISION_LENGTH))->getChildren()->size(),2);
}

TEST_F(TreeTest, testBuildQTree){
    Tree *qTree = tree->BuildQtree();
    tree->printBrokenEdges();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
