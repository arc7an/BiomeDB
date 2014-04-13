angular.module("main", []).
    // draws result graph (first level nodes only)
    controller("svg_controller", function($scope) {
        $scope.circles = generate_xy(10, 500, 300);
        $scope.links = connect_circles($scope.circles);
    }).
    /*
        Model structure. By clicking on the cross to open links and nodes, related to the given node, we make a new request to server.
        We pass data for links and data for nodes.
        $scope.node[k].data = { ... }
        $scope.node[k].neighbours = [{link: lnk1_data}, {node: node1_data}, {link: lnk2_data}, {node: node2_data}]
        where k is an index of the node
     */
    controller("restful_exchange", function($scope, $http) {

        // extract id from "self" string in json data
        var get_id = function(str) {
                return str.split("/").last();
        }

        // returns folded/unfolded icon url depending on the current folding state
        $scope.get_tree_img_src = function(node) {
            var node_id = node.id;
            if($scope.tree[node_id] == true) return "img/minus.gif";
            else return "img/plus.gif";
        }

        // output edges data
        $scope.show_rel = function(node) {
            if(node.rel.type) {
                return node.rel.type + '&mdash';
            }
            return "";
        }

        // swaps the folding state of the given node
        $scope.swap_folding = function(node) {
            var node_id = node.id;
            // TODO: for some reason the node sometimes folds at clicking instead of unfolding
            $scope.tree[node_id] = !$scope.tree[node_id];
        }

        $scope.show_node_data = function(node) {
            $scope.current_node = node;
        }

        $scope.load_children = function(node) {
            var node_id = node.id;
            // if node is not unfolded and its content is not loaded yet
            if(node.nodes.length == 0) {
                $http.post("http://biocastle.net:7474/db/data/cypher", {"query": "start n=node(" + node_id +") match (n)-[r]->(nn) return nn, r"})
                .success(function(data) {
                    console.log(data);
                    for(var i=0; i<data.data.length; i++) {
                        // fill data on vertexes (nodes)
                        node.nodes[i] = {};
                        node.nodes[i].data = data.data[i][0].data;
                        node.nodes[i].nodes = [];
                        var id = get_id(data.data[i][0].self);
                        node.nodes[i].id = id;
                        // fill data on edges
                        node.nodes[i].rel = {};
                        node.nodes[i].rel.type = data.data[i][1].type;
                        // the element is folded by default
                        $scope.tree[id] = false;
                    }
                })
            }
        }

        $scope.send = function() {
            $http.post("http://biocastle.net:7474/db/data/cypher", {"query": $scope.query})
                .success(function(data) {
                    // array of all loaded nodes of first level:
                    $scope.nodes = [];
                    // array to control folded-unfolded state of the node by its id in its interface representation:
                    $scope.tree = [];
                    // using just first column in query response
                    var c = 0;
                    for(var i=0; i<data.data.length; i++) {
                        $scope.nodes[i] = {};
                        $scope.nodes[i].data = data.data[i][c].data;
                        $scope.nodes[i].nodes = [];
                        var id = get_id(data.data[i][0].self);
                        $scope.nodes[i].id = id;
                        $scope.tree[id] = false;
                    }
                })
        }
    })
    // to be able just to press Enter to send cypher query written in the input field
    .directive("ngEnter", function() {
        return function(scope, element, attrs) {
            element.bind("keydown keypress", function(event) {
                if(event.which === 13) {
                    scope.$apply(function(){
                        scope.$eval(attrs.ngEnter)
                    });

                    event.preventDefault();
                }
            });
        };
    })
    .filter("draw_link", function() {
        return function(input) {
            if(input) return '&mdash;' + input + '&mdash;&gt;';
        }
    })
