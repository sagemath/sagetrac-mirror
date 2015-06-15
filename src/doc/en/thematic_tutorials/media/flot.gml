graph [
version 2
Creator "nthiery"
options [
allow_multiple_edges 1
allow_self_loops 1
allow_bends 1
autonumber_nodes 0
autonumber_nodes_by_id 0
autonumber_nodes_by_degree 0
autonumber_edges 0
autonumber_edges_by_id 0
]
directed 1
node_style [
name "default_node_style"
style [
graphics [
w 16.0
h 16.0
]
]
]
edge_style [
name "default_edge_style"
style [
graphics [
]
]
]
node [
id 9
label "s"
graphics [
x 120.0
y 60.0
w 14.0
]
LabelGraphics [
type "text"
]
]
node [
id 12
label "i1"
graphics [
x 260.0
y 60.0
w 18.0
]
LabelGraphics [
type "text"
]
]
node [
id 15
label "i2"
graphics [
x 200.0
y 140.0
w 18.0
]
LabelGraphics [
type "text"
]
]
node [
id 18
label "p"
graphics [
x 400.0
y 60.0
w 15.0
]
LabelGraphics [
type "text"
]
]
node [
id 21
label "i3"
graphics [
x 260.0
y 140.0
w 18.0
]
LabelGraphics [
type "text"
]
]
node [
id 24
label "i4"
graphics [
x 320.0
y 140.0
w 18.0
]
LabelGraphics [
type "text"
]
]
edge [
source 9
target 12
label "&lt;=2"
graphics [
arrow "last"
]
LabelGraphics [
type "text"
]
]
edge [
source 9
target 15
label "&lt;=3"
graphics [
arrow "last"
]
LabelGraphics [
type "text"
]
]
edge [
source 15
target 21
label "&lt;=2"
graphics [
arrow "last"
]
LabelGraphics [
type "text"
]
]
edge [
source 12
target 21
label "&lt;=4"
graphics [
arrow "last"
]
LabelGraphics [
type "text"
]
]
edge [
source 12
target 18
graphics [
arrow "last"
]
]
edge [
source 21
target 24
label "&lt;=1"
graphics [
arrow "last"
]
LabelGraphics [
type "text"
]
]
edge [
source 24
target 18
label "&lt;=7"
graphics [
arrow "last"
]
LabelGraphics [
type "text"
]
]
]

