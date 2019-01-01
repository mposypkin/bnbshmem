/*
 * A simple multithreaded interval-based bnb solver
 */

/*
 * File:   tutorialbnb.cpp
 * Author: yamchenko.y.v
 *
 * Created on January 3, 2018, 5:10 PM
 */

#include <iostream>
#include <limits>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include <functional>
#include <thread>
#include <chrono>
#include <forward_list>
#include <mutex>
#include <memory>
#include <atomic>
#include <unordered_set>
#include <condition_variable>

#include <common/parbench.hpp>
#include <common/bnbstat.hpp>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

constexpr char gKnownRecord[] = "knrec";

static double gEps = 0.1;

static std::string gKnrec = "knrec";

static double gSubsSplitCoeff = 0.5;

static double gStepsSplitCoeff = 0.5;

static int gMtSubsLimit = 10;

static int gMaxStepsTotal = 1000000;

static int gProcs = 4;

static int gMtStepsLimit = 1000;

std::atomic<double> gRecord;

std::vector<BnBStat> stat;

/**
 * BnB state
 */
struct State {

    enum Status {
        IS_READY, IS_PROCESSING, IS_FINISHED
    };

    State(Status stat = IS_READY) : mStatus(stat) {
    }

    State(const State& st) : mRecordVal(st.mRecordVal),
    mRecord(st.mRecord),
    mPool(st.mPool),
    mRemainingSteps(st.mRemainingSteps),
    mStatus(st.mStatus.load()),
    mMutex() {
    }

    State& operator=(const State& st) {
        mRecordVal = st.mRecordVal;
        mRecord = st.mRecord;
        mPool.assign(st.mPool.begin(), st.mPool.end());
        mRemainingSteps = st.mRemainingSteps;
        mStatus.store(st.mStatus.load());
        return *this;
    }

    void merge(const std::shared_ptr<State> s) {
        mRemainingSteps += s->mRemainingSteps;
        s->mRemainingSteps = 0;
        if (s->mRecordVal < mRecordVal) {
            mRecordVal = s->mRecordVal;
            mRecord = s->mRecord;
        }
        mPool.insert(mPool.end(), s->mPool.begin(), s->mPool.end());
    }

    bool hasResources() {
        return !mPool.empty() && (mRemainingSteps > 0);
    }

    void printResourses() {
        std::cout << "Remaining task count: " << mPool.size() << "\n"
                "Remainng steps count: " << mRemainingSteps << "\n";
    }

    void assignTaskTo(std::shared_ptr<State> st) {
        *st = *this;
        mRemainingSteps = 0;
        mPool.clear();
    }

    std::shared_ptr<State> trySplit(std::shared_ptr<State> s = nullptr) {
        std::lock_guard<std::mutex> lock(mMutex);
        const int pool_size = mPool.size();
        if (pool_size <= gMtSubsLimit || mRemainingSteps <= gMtStepsLimit) {
            return nullptr;
        }
        if (s == nullptr) {
            s = std::make_shared<State>();
        }
        const int mv_step_count = mRemainingSteps * gStepsSplitCoeff;
        mRemainingSteps -= mv_step_count;
        s->mRemainingSteps = mv_step_count;
        s->mRecord = mRecord;
        s->mRecordVal = mRecordVal;
        const int mv_sub_count = pool_size * gSubsSplitCoeff;
        auto mv_end_iter = mPool.begin() + mv_sub_count;
        s->mPool.assign(mPool.begin(), mv_end_iter);
        mPool.erase(mPool.begin(), mv_end_iter);
        return s;
    }

    Box getSub() {
        std::lock_guard<std::mutex> lock(mMutex);
        const Box box = mPool.back();
        mPool.pop_back();
        mRemainingSteps--;
        return box;
    }

    void clear() {
        mRecordVal = std::numeric_limits<double>::max();
        mRecord.clear();
        mPool.clear();
        mRemainingSteps = 0;
    }

    void setReady() {
        clear();
        mStatus = IS_READY;
    }

    void setProcessing() {
        mStatus = IS_PROCESSING;
    }

    void setFinished() {
        mStatus = IS_FINISHED;
    }

    bool isFinished() const {
        return (this->mStatus == IS_FINISHED) ? true : false;
    }

    bool isProcessing() const {
        return (this->mStatus == IS_PROCESSING) ? true : false;
    }

    bool isReady() const {
        return (this->mStatus == IS_READY) ? true : false;
    }

    double mRecordVal;
    std::vector<double> mRecord;
    std::vector<Box> mPool;
    int mRemainingSteps;
    std::atomic<Status> mStatus;
    std::mutex mMutex;
};

class Notifier {
public:

    Notifier() : mNotificationCount(0) {
    }

    void notify() {
        mCV.notify_one();
        std::atomic_fetch_add(&mNotificationCount, 1);
    }

    void resolve() {
        int nc = mNotificationCount.load();
        if (nc > 0) {
            std::atomic_fetch_add(&mNotificationCount, 1);
        }
    }

    void wait_notification(std::shared_ptr<State> s) {
        std::unique_lock<std::mutex> lock(s->mMutex);
        mCV.wait_for(lock, mWaitPeriod, [&]() {
            return mNotificationCount.load() != 0;
        });
    }

private:
    std::condition_variable mCV;
    std::atomic<int> mNotificationCount;
    std::chrono::duration<int64_t, std::milli> mWaitPeriod = std::chrono::milliseconds(100);
};

Notifier gNotifier;

struct ThreadList {

    std::shared_ptr<State> getReadyState() {
        if (this->mReadyStates.empty()) {
            return nullptr;
        }

        std::shared_ptr<State> rstate = *(this->mReadyStates.end() - 1);
        this->mReadyStates.pop_back();
        return rstate;
    }

    void addReadyState(std::shared_ptr<State> sptr) {
        this->mReadyStates.push_back(sptr);
    }

    std::forward_list<std::shared_ptr<State>>& operator()() {
        return mList;
    }

    size_t Size() const {
        return mThreadListSize;
    }

    size_t ActiveCount() const {
        return mActiveThreadCount;
    }

    void PushFront(std::shared_ptr<State> st) {
        ++mThreadListSize;
        mList.push_front(st);
    }

    void clear() {
        mList.clear();
        this->mReadyStates.clear();
        mThreadListSize = 0;
        mActiveThreadCount = 0;
    }

    std::forward_list<std::shared_ptr<State>> mList;
    std::vector<std::shared_ptr<State>> mReadyStates;
    size_t mThreadListSize = 0;
    size_t mActiveThreadCount = 0;
};

ThreadList gThreadList;

std::ostream& operator<<(std::ostream & out, const std::shared_ptr<State> s) {
    out << "\"recval\" : " << s->mRecordVal << "\n";
    out << "\"record\" : [";
    for (int i = 0; i < s->mRecord.size(); i++) {
        out << s->mRecord[i];
        if (i != s->mRecord.size() - 1)
            out << ", ";
    }
    out << "]\n";
    out << "\"max steps\" :" << s->mRemainingSteps << "\n";
    return out;
}

double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::shared_ptr<State> s) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return len(f) < len(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = len(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;

    // LOCK AT EVERY ITERATION !!!
    std::lock_guard<std::mutex> lock(s->mMutex);
    s->mPool.push_back(std::move(b1));
    s->mPool.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

void solveSerial(std::shared_ptr<State> s, const BM& bm) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);

    while (s->hasResources()) {
        Box b = s->getSub();
        getCenter(b, c);
        double v = bm.calcFunc(c);
        double rv = gRecord.load();
        while (v < rv) {
            if (gRecord.compare_exchange_strong(rv, v)) {
                s->mRecordVal = v;
                s->mRecord = c;
            }
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= gRecord.load() - gEps) {
            split(b, s);
        }
    }
    s->setFinished();
    gNotifier.notify();
}

void runThread(std::shared_ptr<State> s, const BM& bm) {
    gThreadList.mActiveThreadCount++;
    s->setProcessing();
    std::thread init_thread(solveSerial, s, std::ref(bm));
    init_thread.detach();
}

bool try_assign_run(std::shared_ptr<State> sender, std::shared_ptr<State> receiver, const BM& bm) {
    if (sender->hasResources() && sender->mRemainingSteps > gMtStepsLimit) {
        sender->assignTaskTo(receiver);
        runThread(receiver, bm);
        return true;
    } else {
        return false;
    }
}

void solve(std::shared_ptr<State> init_s, const BM& bm) {
    gRecord.store(std::numeric_limits<double>::max());
    std::shared_ptr<State> first_s = std::make_shared<State>();
    init_s->assignTaskTo(first_s);
    gThreadList.PushFront(first_s);

    runThread(first_s, bm);

    while (gThreadList.ActiveCount() || (init_s->hasResources() && init_s->mRemainingSteps > gMtStepsLimit)) {
        for (auto iter = gThreadList().begin(); iter != gThreadList().end(); ++iter) {
            std::shared_ptr<State> cur_state = *iter;
            if (cur_state->isReady()) {
                continue;
            }

            if (cur_state->isFinished()) {
                init_s->merge(cur_state);
                cur_state->setReady();
                gThreadList.mActiveThreadCount--;
                gNotifier.resolve();

                bool flag = try_assign_run(init_s, cur_state, bm);
                if (!flag) {
                    gThreadList.addReadyState(cur_state);
                }
                continue;
            }

            if (gThreadList.mActiveThreadCount < gProcs) {
                std::shared_ptr<State> new_state;
                if (gThreadList.mActiveThreadCount == gThreadList.Size()) {
                    new_state = cur_state->trySplit();
                    if (new_state == nullptr)
                        continue;

                    gThreadList.PushFront(new_state);
                } else {
                    std::shared_ptr<State> ready_state = gThreadList.getReadyState();
                    if (ready_state == nullptr) {
                        std::cout << "Can't find ready. Ready size : " << gThreadList.mReadyStates.size() << std::endl;
                        continue;
                    }

                    new_state = cur_state->trySplit(ready_state);
                    if (new_state == nullptr) {
                        gThreadList.addReadyState(ready_state);
                        continue;
                    }
                }
                runThread(new_state, bm);
            } else {
                gNotifier.wait_notification(init_s);
            }
        }
    }

    std::cout << "Remaining steps: " << init_s->mRemainingSteps << std::endl;
    if (init_s->hasResources()) {
        init_s->printResourses();
        solveSerial(init_s, bm);
    }

    gThreadList.clear();

    std::cout << "Problem is solved" << std::endl;
}

double findMin(const BM& bm) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    std::shared_ptr<State> s = std::make_shared<State>();
    s->mPool.push_back(ibox);
    s->mRecordVal = std::numeric_limits<double>::max();
    s->mRemainingSteps = gMaxStepsTotal;

    if (gKnrec == gKnownRecord) {
        s->mRecordVal = bm.getGlobMinY();
    } else {
        s->mRecordVal = std::numeric_limits<double>::max();
    }

    std::chrono::time_point<std::chrono::steady_clock> start, end;
    start = std::chrono::steady_clock::now();
#if 0
    solveSerial(&s, bm);
#else
    solve(s, bm);
#endif
    end = std::chrono::steady_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    const int steps = gMaxStepsTotal - s->mRemainingSteps;
    std::cout << "Time per subproblem: " << ((mseconds > 0) ? ((double) mseconds / (double) steps) : 0) << " miscroseconds." << std::endl;
    if (s->mRemainingSteps == 0) {
        std::cout << "Failed to converge in " << gMaxStepsTotal << " steps\n";
    } else {
        std::cout << "Converged in " << steps << " steps\n";
    }

    std::cout << "BnB found = " << gRecord.load() << std::endl;
    std::cout << " at x [ ";
    std::copy(s->mRecord.begin(), s->mRecord.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";

    stat.emplace_back((double) mseconds, steps);

    return gRecord.load();
}

void printInfo() {
    std::cout << "Record is lock-free: " << gRecord.is_lock_free() << std::endl;
    if (gKnrec == gKnownRecord) {
        std::cout << "Init global record: true" << std::endl;
    } else {
        std::cout << "Init global record: false" << std::endl;
    }
    std::cout << "Eps: " << gEps << std::endl;
    std::cout << "Max steps count: " << gMaxStepsTotal << std::endl;
    std::cout << "Max thread count: " << gProcs << std::endl;
    std::cout << "Available thread count: " << std::thread::hardware_concurrency() << std::endl;
    std::cout << "Steps split coeffitient count: " << gStepsSplitCoeff << std::endl;
    std::cout << "Subs split coeffitient count: " << gSubsSplitCoeff << std::endl;
    std::cout << "Steps limit: " << gMtStepsLimit << std::endl;
    std::cout << "Subs limit: " << gMtSubsLimit << std::endl << std::endl;
}

bool testBench(const BM& bm) {
    std::cout << "*************Testing benchmark**********" << std::endl;
    printInfo();
    std::cout << bm;

    bool res_flag = true;
    double v = findMin(bm);
    double diff = v - bm.getGlobMinY();
    if (diff > gEps) {
        res_flag = false;
        std::cout << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;

    return res_flag;
}



int main(int argc, char** argv) {
    std::string bench;
#if 0    
    Benchmarks<double> tests;
#else 
    ParBenchmarks<double> tests;
#endif    
    int nruns = 0;

    if ((argc == 2) && (std::string(argv[1]) == std::string("list"))) {
        for (auto b : tests) {
            std::cout << b->getDesc() << "\n";
        }
        return 0;
    } else if (argc == 11) {
        nruns = atoi(argv[1]);
        bench = argv[2];
        gKnrec = argv[3];
        gEps = atof(argv[4]);
        gMaxStepsTotal = atoi(argv[5]);
        gProcs = atoi(argv[6]);
        gStepsSplitCoeff = atof(argv[7]);
        gSubsSplitCoeff = atof(argv[8]);
        gMtStepsLimit = atoi(argv[9]);
        gMtSubsLimit = atoi(argv[10]);
    } else {
        std::cerr << "Usage: " << argv[0] << "num_of_runs name_of_bench knrec|unknrec eps max_steps thread_count steps_split_coeff subs_split_coeff step_limit subs_limit\n";
        std::cerr << "or to list benchmarks run:\n";
        std::cerr << argv[0] << " list\n";
        return -1;
    }
    for (int z = 0; z < nruns; z++)
        for (auto bm : tests) {
            if (bench == bm->getDesc())
                testBench(*bm);
        }
    std::cout << "Statistics:\n" << stat;
}
