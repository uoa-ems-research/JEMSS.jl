# common definitions

const Float = Float64
# const Int = Int64
# const UInt = UInt64

# float and int types for storing all-pairs shortest paths data
# can reduce precision in order to reduce memory used
const FloatSpTime = Float64 # precision for storing shortest path travel times
const IntRNode = Int16 # precision for storing number of nodes in rGraph
const IntFadj = Int8 # precision for storing maximum number of nodes adjacent to any node in rGraph (= length of longest vector in network.rGraph.fadjList)

const sourcePath = splitdir(@__FILE__)[1]

# run modes
const debugMode = false
const checkMode = true # for data checking, e.g. assertions that are checked frequently

# file chars
const delimiter = ';'
const newline = "\r\n"

# misc null values
const nullIndex = -1
const nullX = -200 # outside lat/lon range
const nullY = -200
const nullTime = -1.0
const nullDist = -1.0

# call/travel priorities
@enum Priority nullPriority=0 highPriority=1 medPriority=2 lowPriority=3

# ambulance classes
@enum AmbClass nullAmbClass=0 als=1 bls=2 # als = advanced life support, bls = basic life support

# event types/forms
@enum EventForm nullEvent ambGoesToSleep ambWakesUp callArrives considerDispatch ambDispatched ambReachesCall ambGoesToHospital ambReachesHospital ambBecomesIdle ambReachesStation ambRedirected considerMoveUp ambMoveUp

# ambulance statuses
@enum AmbStatus ambNullStatus ambSleeping ambIdleAtStation ambGoingToCall ambAtCall ambGoingToHospital ambAtHospital ambGoingToStation ambMovingUp

# call statuses
@enum CallStatus callNullStatus callScreening callQueued callWaitingForAmb callOnSceneCare callGoingToHospital callAtHospital callProcessed

@enum MoveUpModule nullMoveUpModule compTableModule dmexclpModule priorityListModule zhangIpModule temp1Module temp2Module
