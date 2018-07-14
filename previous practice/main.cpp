#include <iostream>
#include <random>
#include <time.h>
#include <iomanip>
#include <vector>
#include <deque>
#include <queue>
#include <algorithm>
#include <fstream>
#include <string>


// train system status
#define ARRIVAL 0
#define WAITING 1
#define DOCKED 2
#define HOGGED_DS 3
#define HOGGED_IQ 4
#define CREW_ARRIVES 5
#define DONE 6
#define KILL -1

// default parameters
#define DEFAULT_SIMULATION_TIME 72000
#define DEFAULT_INTERVAL 10
#define DEFAULT_CREW_SHIFT 12


// =====================================================================
// STRUCTURES 
// =====================================================================

struct Dock
{
	Dock() : isDockOpen(true), dockBusyTime(0.0), dockHogoutTime(0.0)  {}

	bool isDockOpen;
	double dockBusyTime;
	double dockHogoutTime;
};

struct Train
{
	Train(double at, double ut, double ch)
		: arrivalTime(at), unloadingTime(ut), crewHours(ch), 
		nextEventTime(at), previousEventTime(at), timeInQueue(0.0), 
		timeInSystem(0.0), hogoutCount(0), status(ARRIVAL), docked(false) {}

	int trainID;
	int crewID;

	double arrivalTime;
	double unloadingTime;
	double crewHours;
	double nextEventTime; 
	double previousEventTime;
	double timeInQueue;
	double timeInSystem; // time spent from arrival to departure
	int hogoutCount;
	int status; // used to determine train state in the FSM-ish logic 
	bool docked; // used to allow a crew into the docks if previous crew hogged
};

// =====================================================================
// SIMULATION MESSAGES
// =====================================================================

void arrivalMessage(const Train& t, int queueNum)
{
	std::cout << std::fixed << "Time " << t.nextEventTime
		<< ": train " << t.trainID
		<< " arrival for " << t.unloadingTime << "h of unloading, " << std::endl;

	std::cout << "\t\tcrew " << t.crewID
		<< " with " << t.crewHours << "h before hogout "
		<< "(Q=" << queueNum << ")" << std::endl << std::endl;
}

void departureMessage(const Train& t)
{
	std::cout << std::fixed << "Time " << t.nextEventTime
		<< ": train " << t.trainID << " departing" << std::endl << std::endl;
}

void dockingMessage(const Train& t)
{
	std::cout << std::fixed << "Time " << t.nextEventTime
		<< ": train " << t.trainID
		<< " entering dock for " << t.unloadingTime << "h of unloading, " << std::endl;

	std::cout << "\t\tcrew " << t.crewID
		<< " with " << t.crewHours << "h before hogout " << std::endl << std::endl;
}

void hogoutMessageIQ(const Train& t)
{
	std::cout << std::fixed << "Time " << t.nextEventTime 
		<< ": train " << t.trainID << " crew " << t.crewID
		<< " hogged out in queue (SERVER HOGGED)" << std::endl << std::endl;
}

void hogoutMessageDS(const Train& t)
{
	std::cout << std::fixed << "Time " << t.nextEventTime 
		<< ": train " << t.trainID << " crew " << t.crewID 
		<< " hogged out during service (SERVER HOGGED) " << std::endl << std::endl;
}

void newCrewMessage(const Train& t)
{
	std::cout << std::fixed << "Time " << t.nextEventTime 
		<< ": train " << t.trainID
		<< " replacement crew " << t.crewID << " arrives (SERVER UNHOGGED)" << std::endl << std::endl;
}

void displayStats(const std::vector<int>& histogram, const Dock& dock, 
	int numTrains, int maxQueue, double avgTIS, 
	double maxTIS, double avgTIQ, double elapsedTime)
{
	std::cout << std::endl << std::endl << std::endl
		<< "Statistics" << std::endl
		<< "-------------------------" << std::endl;

	std::cout << "Total number of trains served: "
		<< numTrains << std::endl;

	std::cout << "Average time-in-system per train: "
		<< avgTIS / numTrains << "h" << std::endl;

	std::cout << "Maximum time-in-system per train: "
		<< maxTIS << "h" << std::endl;

	std::cout << "Average time-in-queue over trains: "
		<< avgTIQ / numTrains << "h" << std::endl;

	std::cout << "Dock busy percentage: " << std::setprecision(2)
		<< 100 - (dock.dockBusyTime / elapsedTime) * 100 << "%" << std::endl;

	std::cout << "Dock idle percentage: " << std::setprecision(2)
		<< (dock.dockBusyTime / elapsedTime) * 100 << "%" << std::endl;

	std::cout << "Dock hogged-out percentage: " << std::setprecision(2)
		<< (dock.dockHogoutTime / elapsedTime) * 100 << "%" << std::endl;

	std::cout << "Maximum number of trains in queue: "
		<< maxQueue << std::endl;

	std::cout << "Histogram of hogout count per train:" << std::endl;

	for (int i = 0; i < histogram.size(); ++i)
	{
		if (histogram[i] == 0) continue;
		else std::cout << "[" << i << "]: " << histogram[i] << std::endl;
	}
}

void endSimulationMessage(double elapsedTime)
{
	std::cout << std::fixed << "Time " << elapsedTime << ": "
		<< "simulation ended" << std::endl;
}

// =====================================================================
// SIMULATION HANDLERS
// =====================================================================

void handleArrival(Train& t, int& numTrains, int& numCrews,
	int& queueNum, int& maxQueue, double simulationTime)
{
	t.status = WAITING;
	if (++queueNum > maxQueue) maxQueue = queueNum;
	t.trainID = numTrains++;
	t.crewID = numCrews++;
	t.timeInSystem = t.arrivalTime;
	arrivalMessage(t, queueNum);
}

void handleWaiting(Train& t, int& prevTrain, double& time, Dock& dock)
{
	// only lets train in dock if 1) the dock is open 
	// and 2) the train is next in line, compares train IDs
	if (dock.isDockOpen && t.trainID - 1 == prevTrain)
	{
		t.status = DOCKED;
		dock.isDockOpen = !dock.isDockOpen;
	}
	else //subtract crew hours while waiting for dock to open
	{
		// the 0.1 is to increment the event time in a case where there is no
		// available or appropriate time to calculate the next event
		if (time < t.previousEventTime) time = t.previousEventTime + 0.1;
		double temp = time - t.previousEventTime;
		if (temp == 0) temp = 0.1;

		if (t.crewHours - temp < 0)
		{
			t.status = HOGGED_IQ;
			++t.hogoutCount;
			t.previousEventTime = t.nextEventTime;
			t.nextEventTime += t.crewHours;
			t.timeInQueue += t.crewHours;
			time = t.nextEventTime;
		}
		else
		{
			t.crewHours -= temp;
			t.previousEventTime += temp;
			t.timeInQueue += temp;
			t.nextEventTime = t.previousEventTime;
		}
	}
}

void handleDocked(Train& t, Dock& dock, double& time)
{
	if (!t.docked)
	{
		t.docked = !t.docked;
		dockingMessage(t);
	}

	if (t.crewHours >= t.unloadingTime)
	{
		t.status = DONE;
		t.previousEventTime = t.nextEventTime;
		t.nextEventTime += t.unloadingTime;
		dock.dockBusyTime += t.unloadingTime;
		time = t.nextEventTime;
	}
	else // if crew hours left < required unloading time 
	{
		t.status = HOGGED_DS;
		++t.hogoutCount;
		t.unloadingTime -= t.crewHours;
		t.previousEventTime = t.nextEventTime;
		t.nextEventTime += t.crewHours;
		dock.dockBusyTime += t.crewHours;
		time = t.nextEventTime;
	}
}

void handleHoggedIQ(Train& t, int& numCrews, 
	double& time, std::queue<double> travelTimes)
{
	t.status = CREW_ARRIVES;
	double travTime = travelTimes.front();
	travelTimes.pop();
	t.crewHours = 12 - travTime;
	t.crewID = numCrews++;
	t.previousEventTime = t.nextEventTime;
	t.nextEventTime += travTime;
	t.timeInQueue += travTime;
	time = t.nextEventTime;
	hogoutMessageIQ(t);
}

void handleHoggedDS(Train& t, Dock& dock, int& numCrews, 
	double& time, std::queue<double> travelTimes)
{
	t.status = CREW_ARRIVES;
	double travTime = travelTimes.front();
	travelTimes.pop();
	t.crewHours = 12 - travTime;
	t.crewID = numCrews++;
	t.previousEventTime = t.nextEventTime;
	t.nextEventTime += travTime;
	dock.dockHogoutTime += travTime;
	time = t.nextEventTime;
	hogoutMessageDS(t);
}

void handleCrewArrives(Train& t)
{
	t.status = (t.docked) ? DOCKED : WAITING;
	newCrewMessage(t);
}

void handleDone(Train& t, int& prevTrain, int& numQueue, Dock& dock)
{
	t.status = KILL; // kill status, used in the main loop to pop off queue
	departureMessage(t);
	--numQueue;
	++prevTrain;
	t.docked = !t.docked;
	dock.isDockOpen = !dock.isDockOpen;
	t.timeInSystem = t.nextEventTime - t.arrivalTime;
}

void run(Train& t, int& numTrains, int& numCrews, 
	int& numQueue, int& maxQueue, int& prevTrain, double& time, 
	double simulationTime, Dock& dock, std::queue<double> travelTimes)
{
	// handled train states like a finite state machine with 6 possible states
	switch (t.status)
	{
	case ARRIVAL:
		handleArrival(t, numTrains, numCrews, 
			numQueue, maxQueue, simulationTime); break;
	case WAITING:
		handleWaiting(t, prevTrain, time, dock); break;
	case DOCKED:
		handleDocked(t, dock, time); break;
	case HOGGED_IQ:
		handleHoggedIQ(t, numCrews, time, travelTimes); break;
	case HOGGED_DS:
		handleHoggedDS(t, dock, numCrews, time, travelTimes); break;
	case CREW_ARRIVES:
		handleCrewArrives(t); break;
	case DONE:
		handleDone(t, prevTrain, numQueue, dock); break;
	}
}

// =====================================================================
// HELPER FUNCTIONS
// =====================================================================

// comparator function used to sort the deque by event times
// if times are equal, sort by lowest train ID
bool comparator(const Train& t1, const Train& t2)
{
	if (t1.nextEventTime == t2.nextEventTime)
		return t1.trainID < t2.trainID;
	else return t1.nextEventTime < t2.nextEventTime;
}

// fill the train arrival schedule if not predetermined
void produceTimes(std::queue<double>& travelTimes, std::mt19937 generator)
{
	std::uniform_real_distribution<double> replacementDist(2.5, 3.5);
	for (int i = 0; i < 500; ++i)
		travelTimes.push(replacementDist(generator));
}

void initPredetermined(std::deque<Train>& trains, std::queue<double>& travelTimes, 
	double& simulationTime, std::string s, std::string t)
{
	std::ifstream infile;
	s = "schedule.txt";
	t = "traveltimes.txt";
	infile.open(s);
	while (!infile.eof())
	{
		double arrival;
		double unloading;
		double hours;
		infile >> arrival >> unloading >> hours;
		trains.push_back(Train(arrival, unloading, hours));
	}
	simulationTime = trains.back().arrivalTime * trains.back().arrivalTime;
	infile.close();
	infile.open(t);
	while (!infile.eof())
	{
		double trav;
		infile >> trav;
		travelTimes.push(trav);
	}
	infile.close();
}

void initRandom(std::deque<Train>& trains, std::mt19937 generator,
	double interarrivalTimes, double simulationTime)
{
	std::exponential_distribution<double> arrivalDist(1.0 / interarrivalTimes);
	std::uniform_real_distribution<double> unloadDist(3.5, 4.5);
	std::uniform_real_distribution<double> crewDist(6.0, 11.0);

	for (double i = 0; i < simulationTime;)
	{
		double arrival = arrivalDist(generator);
		double unloading = unloadDist(generator);
		double hours = crewDist(generator);
		trains.push_back(Train(i += arrival, unloading, hours));
	}
}

void statistics(std::vector<Train>& stats, const Dock& dock,
	int numTrains, int maxQueue, double elapsedTime)
{
	std::vector<int> histogram;
	stats.shrink_to_fit();

	double avgTIQ = 0.0; // avg time in queue = time in queue / #trains
	double avgTIS = 0.0; // avg time in sys = (time in queue + docked + hogged) / #trains
	double maxTIS = 0.0;

	// arbitrary overkill vector size
	for (int i = 0; i < 500; ++i) histogram.push_back(0);

	for (int i = 0; i < stats.size(); ++i)
	{
		++histogram[stats[i].hogoutCount];
		avgTIQ += stats[i].timeInQueue;
		avgTIS += stats[i].timeInSystem;
		if (stats[i].timeInSystem > maxTIS) maxTIS = stats[i].timeInSystem;
	}

	displayStats(histogram, dock, numTrains, maxQueue,
		avgTIS, maxTIS, avgTIQ, elapsedTime);
}

// =====================================================================
// MAIN FUNCTION
// =====================================================================

int main(int argc, char* argv[])
{
	double interarrivalTimes = DEFAULT_INTERVAL;
	double simulationTime = DEFAULT_SIMULATION_TIME;

	if (argc == 3)
	{
		interarrivalTimes = atoi(argv[1]);
		simulationTime = atoi(argv[2]);
	}

	std::mt19937 generator(time(NULL));
	std::deque<Train> trains; // used to store the list of generated trains 
	std::queue<double> travelTimes;

	(argc == 4)
		? initPredetermined(trains, travelTimes, simulationTime, argv[2], argv[3])
		: initRandom(trains, generator, interarrivalTimes, simulationTime);


	std::deque<Train> tq; // used for the train queue of the simulation
							// I chose a deque because you can use subscript indexing
	std::vector<Train> stats; // used to gather the analytics data at the very end;
	Dock dock;
	int numTrains = 0;
	int numCrews = 0;
	int numQueue = -1; 
	int maxQueue = 0;
	int prevTrain = -1;
	double time = trains[0].arrivalTime;
	double elapsedTime = 0.0;

	// MAIN LOOP
	tq.push_back(trains[0]);
	trains.pop_front();

	do
	{
		if (travelTimes.empty()) produceTimes(travelTimes, generator);

		// add new train(s) to queue 
		for (int i = 0; i < trains.size(); ++i)
		{
			if (tq[0].nextEventTime > trains[i].arrivalTime)
			{
				tq.push_back(trains[i]);
				trains.pop_front();
			}
			else
			{
				tq.push_back(trains[0]);
				trains.pop_front();
				break;
			}
		}

		sort(tq.begin(), tq.end(), comparator);
		run(tq[0], numTrains, numCrews, numQueue, maxQueue,
			prevTrain, time, simulationTime, dock, travelTimes);

		if (time > elapsedTime) elapsedTime = time;
		if (tq[0].status == KILL)
		{
			// prevent a train from arriving after the time limit if possible
			if (tq[1].nextEventTime > simulationTime && tq[1].status != DONE) break;
			stats.push_back(tq[0]);
			tq.pop_front();
		}
		if (tq.empty()) break;
	} while (elapsedTime < simulationTime);


	// finish unloading trains that have made it before the time limit
	while (!tq.empty())
	{
		sort(tq.begin(), tq.end(), comparator);

		if (tq[0].status == ARRIVAL) // another check for late arrivals just in case
		{
			tq.pop_front();
			continue;
		}

		run(tq[0], numTrains, numCrews, numQueue, maxQueue,
			prevTrain, time, simulationTime, dock, travelTimes);

		if (tq[0].status == KILL)
		{
			stats.push_back(tq[0]);
			elapsedTime = tq[0].nextEventTime;
			tq.pop_front();
		}
	}

	endSimulationMessage(elapsedTime);
	statistics(stats, dock, numTrains, maxQueue, elapsedTime);

	return 0;
}



// =====================================================================
// TESTING AND DEBUGGING 
// =====================================================================

// some debugging code to display values at each step
//std::cout << "t" << t.trainID << std::endl
//	<< " hours " << t.crewHours << std::endl
//	<< " prev " << t.previousEventTime << std::endl
//	<< " next " << t.nextEventTime << std::endl
//	<< " time " << time << std::endl
//	<< std::endl << std::endl;


// used to iterate through each timestep one by one 
// so I can manually check the values for correctness
//std::string s;
//std::getline (std::cin, s);
