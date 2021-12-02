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
            subTreeRoots->push_back(root);

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
            subTreeRoots->push_back(tree->getTask(2+(int)((DIVISION_LENGTH+1)/3)));
            subTreeRoots->push_back(tree->getTask(2+(int)(2*(DIVISION_LENGTH+1)/3)));

            // build chain part
            // one edge by arount 80% of the chain will be broken
            Task * parent;
            parent = tree->getTask(2);
            subTreeRoots->push_back(parent);
            task = tree->getTask(curr_task_id);
            task->setParent(parent);
            parent->addChild(task);  
            curr_task_id ++;

            
            for (curr_task_id; curr_task_id<(2*DIVISION_LENGTH)+3;curr_task_id++){
                task = tree->getTask(curr_task_id);
                parent =tree->getTask(curr_task_id-1); 
                task->setParent(parent);
                parent->addChild(task); 

                if (curr_task_id == curr_task_id+(int)(((2*DIVISION_LENGTH)+3 - curr_task_id)*0.8)) {
                    task->breakEdge();
                    subTreeRoots->push_back(task);
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

            subTreeRoots->push_back(tree->getTask(curr_task_id));

            for (curr_task_id; curr_task_id<=NUM_TASKS; curr_task_id++){
                parent = tree->getTask((int)((curr_task_id-id_before_binary_tree)/2)+id_before_binary_tree);
                task = tree->getTask(curr_task_id);
                task->setParent(parent);
                parent->addChild(task);  

                if (curr_task_id%4 == 0){
                    task->breakEdge();
                }
                if (curr_task_id%8 == 0){
                   subTreeRoots->push_back(task);
                }
            }
        }

        ~TreeTest() override {
            delete tree;
        }

        bool isChild(vector<Task *> children, int child_id){
            for (Task * cur_child :children){
                if (child_id == cur_child->getId()) return true;
            }
            return false;
        }

    Tree * tree = new Tree();
    vector <Task *> * subTreeRoots = new vector<Task *>();

    };

}

TEST_F(TreeTest, emptyTreeIsEmpty){
    Tree * t = new Tree();
    EXPECT_EQ(t->getSize(),0);
}

// Disclaimer: Only works, when DIVISION_LENGTH is large enough
TEST_F(TreeTest,testFixtureStructure){
    cout << tree->getSize()<<endl;
    EXPECT_EQ(tree->getSize(),4*DIVISION_LENGTH +1);
    EXPECT_EQ(tree->getRoot()->getChildren()->size(),DIVISION_LENGTH +2);
    cout << tree->getTask((int)(1.5*DIVISION_LENGTH))->getId()<<endl;
    EXPECT_EQ(tree->getTask((int)(1.5*DIVISION_LENGTH))->getChildren()->size(),1);
    EXPECT_EQ(tree->getTask((int)(2.5*DIVISION_LENGTH))->getChildren()->size(),2);
}

TEST_F(TreeTest, testBuildQTreeSize){
    Tree *qTree = tree->BuildQtree();

    // there is one node in the QTree for every broken edge in the original tree
    EXPECT_EQ(tree->countBrokenEdges(), qTree->getSize());
}

TEST_F(TreeTest, testBuildQTreeOtherSideIds){
    Tree *qTree = tree->BuildQtree();
    // the root of a broken part has the otherSideId of the next node in the qtree and vice versa
    Task *task;
    Task *qTree_task;
    unsigned long expected_qTree_id = 1;
    for (unsigned long task_id = 1; task_id<tree->getSize()+1; task_id++){
        task = tree-> getTask(task_id);
        

        if (task->isBroken()){
            qTree_task = qTree->getTask(expected_qTree_id);
            EXPECT_EQ(task->getOtherSideId(), expected_qTree_id);
            EXPECT_EQ(task_id, qTree_task->getOtherSideId());
            expected_qTree_id ++;
        }
    }

}

TEST_F(TreeTest, testBuildQTreeRootQuality){
    Tree *qTree = tree->BuildQtree();
    // the root of a qTree has the root attribute
    Task *qTreeRoot = qTree->getTask(tree->getRoot()->getOtherSideId());
    EXPECT_EQ(qTreeRoot->isRoot(), true);
    EXPECT_EQ(qTreeRoot, qTree->getRoot());

}

TEST_F(TreeTest, testBuildQTreeParentRelations){
    Tree *qTree = tree->BuildQtree();
    Task *qParent, *originalParent, *originalTask;
    for (Task *qTask: *qTree->getTasks()){
        if(!qTask->isRoot()){
            qParent = qTask->getParent();
            originalTask = tree->getTask(qTask->getOtherSideId());
            originalParent = tree->getTask(qParent->getOtherSideId());
            EXPECT_TRUE(originalParent->isBroken());
            EXPECT_GT(originalTask->getId(),originalParent->getId());
            EXPECT_TRUE(isChild(*(qParent->getChildren()),qTask->getId()));
        }
    }
}

TEST_F(TreeTest,testBuildQTreeChildRelations){
    Tree *qTree = tree->BuildQtree();
    Task *qChild, *originalChild, *originalTask;
    for (Task *qTask: *qTree->getTasks()){
        for (Task *qChild:*qTask->getChildren()){
            originalTask = tree->getTask(qTask->getOtherSideId());
            originalChild = tree->getTask(qChild->getOtherSideId());
            EXPECT_TRUE(originalChild->isBroken());
            EXPECT_LT(originalTask->getId(),originalChild->getId());
            EXPECT_TRUE(isChild(*(qTask->getChildren()),qChild->getId()));
        }
    }
}

TEST_F(TreeTest, testBuildSubTreeRootAttributesDoNotChange){
    Tree *subTree;
    Task *task;
    Task *subTreeRoot;
    for (int i = 0; i< subTreeRoots->size(); i++){
        task = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, task);
        subTreeRoot = subTree->getRoot();

        EXPECT_EQ(task->getCost(),subTreeRoot->getCost());
        EXPECT_EQ(task->getEdgeWeight(),subTreeRoot->getEdgeWeight());
        EXPECT_EQ(task->getNodeWeight(),subTreeRoot->getNodeWeight());
        EXPECT_EQ(task->getMakespanWeight(),subTreeRoot->getMakespanWeight());
    } 
}

TEST_F(TreeTest, testBuildSubTreeRoots){
    Tree *subTree;
    Task *task;
    Task *subTreeRoot;
    for (int i = 0; i< subTreeRoots->size(); i++){
        task = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, task);
        subTreeRoot = subTree->getRoot();

        EXPECT_EQ(subTreeRoot->isRoot(),true);
        EXPECT_EQ(subTreeRoot->getOtherSideId(),task->getId());
        EXPECT_EQ(task->getOtherSideId(), subTreeRoot->getId());
    }    
}

TEST_F(TreeTest, testBuildSubtreeOtherSideIdFromSubTreeToOriginal){
    Tree *subTree;
    Task *originalTask;
    for (int i = 0; i< subTreeRoots->size(); i++){
        originalTask = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, originalTask);
        for(Task * subTreeTask: *subTree->getTasks()){
            originalTask = tree->getTask(subTreeTask->getOtherSideId());   
            EXPECT_EQ(originalTask->getOtherSideId(), subTreeTask->getId());
        }
    }  
}

TEST_F(TreeTest, testBuildSubtreeOtherSideIdFromOriginalToSubTree){
    Tree *subTree;
    Task *originalTask;
    vector <Task *> *stack ;
    Task * subTreeTask;

    for (int i = 0; i< subTreeRoots->size(); i++){
        originalTask = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, originalTask);
        stack = new vector <Task *>;
        stack->push_back(originalTask);
        while(!stack->empty()){
            originalTask = stack->back();
            stack ->pop_back();
            subTreeTask = subTree->getTask(originalTask->getOtherSideId());
            EXPECT_EQ(subTreeTask->getOtherSideId(), originalTask->getId());
            for(Task *child: *originalTask->getChildren()){
                if(!child->isBroken()){
                    stack->push_back(child);
                }
            }
        }
        
    }  
}

TEST_F(TreeTest, testBuildSubtreeParentRelationships){
    Tree *subTree;
    Task *originalTask;
    Task *originalParent;
    Task *subTreeParent;

     for (int i = 0; i< subTreeRoots->size(); i++){
        originalTask = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, originalTask);
        
        for (Task *subTreeTask : *subTree->getTasks()){
            originalTask = tree->getTask(subTreeTask ->getOtherSideId());
            if (!subTreeTask->isRoot()){
                originalParent = originalTask->getParent();
                subTreeParent = subTreeTask->getParent();
                EXPECT_EQ(originalParent->getOtherSideId(), subTreeParent->getId());
            }
            else{
                EXPECT_EQ(subTreeTask->getParentId(),0);
            }
        }
    }
}

TEST_F(TreeTest, testBuildSubtreeChildRelationships){
    Tree *subTree;
    Task *originalTask;
    Task *subTreeChild;
     for (int i = 0; i< subTreeRoots->size(); i++){
        originalTask = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, originalTask);
        
        for (Task *subTreeTask : *subTree->getTasks()){
            originalTask = tree->getTask(subTreeTask ->getOtherSideId());
            for (Task *originalChild: *originalTask->getChildren()){
                if (!originalChild->isBroken()){
                    int subTreeChildId = originalChild->getOtherSideId();
                    EXPECT_EQ(isChild(*(subTreeTask->getChildren()),subTreeChildId),true);
                }
            }        
        }
    }
}

TEST_F(TreeTest, testBuildSubTreeOthersideIdOfBrokenEdgesIsZero){
    Tree *subTree;
    Task *originalTask;
     for (int i = 0; i< subTreeRoots->size(); i++){
        originalTask = subTreeRoots->at(i);
        subTree = BuildSubtree(tree, originalTask);
        
        for (Task *subTreeTask : *subTree->getTasks()){
            originalTask = tree->getTask(subTreeTask ->getOtherSideId());
            for (Task *originalChild: *originalTask->getChildren()){
                if (originalChild->isBroken()){
                    EXPECT_EQ(originalChild->getOtherSideId(),0);
                }
            }        
        }
    }
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
