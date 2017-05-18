var chartGenerator = {
  Date_TU_Global_China : function  (data) {
    var globalData = [];
    var chinaData = [];

    for(var i = 0 ; i < data["Date_TU_Global"].length; i ++) {
      globalData.push([Date.parse(data["Date_TU_Global"][i].date), parseInt(data["Date_TU_Global"][i]["Total Users"])])
    }

    for(var i = 0 ; i < data["Date_TU_China"].length; i ++) {
      chinaData.push([Date.parse(data["Date_TU_China"][i].date), parseInt(data["Date_TU_China"][i]["Total Users"])])
    }
    $('#Date_TU_Global_China').highcharts({
      colors: ["#0000ff", "#ff0066", "#2b908f", "#90ee7e", "#eeaaee",
    "#55BF3B", "#DF5353", "#7798BF", "#aaeeee"],
        chart: {
      backgroundColor: null,
            zoomType: 'x'
        },
        title: {
            text: 'Total User Growth'
        },
        subtitle: {
            text: 'Global & China'
        },

        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },

              title: {
                  text: 'Date'
              },
              floor: Date.UTC(2015, 00, 01)
          },
          yAxis: {
              title: {
                  text: 'Users'
              },
              min: 0
          },
          tooltip: {
              headerFormat: '<b>{series.name}</b><br>',
              pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
          },
        legend: {
            enabled: true
        },
        plotOptions: {
            area: {
                color: '#0038ff',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            name: 'Global',
            data: globalData
        },
        {
            name: 'China',
            data: chinaData
        }]
    });
  },
  Date_Global_Country_TU : function(data) {
    var rankRegion = [];
    var rankRegionName = [];
    var rankRegionBasic = {}
    rankRegionBasic = {"Others" : 0}
    var others = {}; // sum of other regions,
    var top = {};
    for(var i = 0 ; i < data["Global_Country_TU"].length && i < 4; i ++) { //find top 4
      rankRegion.push(data["Global_Country_TU"][i]);
      rankRegionName.push(data["Global_Country_TU"][i].Country);
      rankRegionBasic[data["Global_Country_TU"][i].Country] =  0;
    }
    for(var i = 0 ; i < data["Date_Global_Country_TU"].length; i++) {
      if(rankRegionName.indexOf(data["Date_Global_Country_TU"][i].country) >= 0 ) { // if this region in top 4
        if(data["Date_Global_Country_TU"][i].country == "TW")
          console.info(data["Date_Global_Country_TU"][i])
        var sum = rankRegionBasic[data["Date_Global_Country_TU"][i].country] + parseInt(data["Date_Global_Country_TU"][i]["Total Users"])
        rankRegionBasic[data["Date_Global_Country_TU"][i].country] = sum
        var arr = [Date.parse(data["Date_Global_Country_TU"][i].date), sum];
        if(top[data["Date_Global_Country_TU"][i].country] && top[data["Date_Global_Country_TU"][i].country].length > 0) {
          top[data["Date_Global_Country_TU"][i].country].push(arr);
        } else {
          top[data["Date_Global_Country_TU"][i].country] = [];
          top[data["Date_Global_Country_TU"][i].country].push(arr)
        }
      } else { //others
        var sum = rankRegionBasic["Others"] + parseInt(data["Date_Global_Country_TU"][i]["Total Users"])
        rankRegionBasic["Others"] = sum
        var arr = [Date.parse(data["Date_Global_Country_TU"][i].date), sum];
        if(top["Others"] && top["Others"].length > 0) {
          top["Others"].push(arr);
        } else {
          top["Others"] = [];
          top["Others"].push(arr)
        }
      }
    }
    var seriesData = []; rankRegionName.push("Others")
    for(var k in rankRegionName) {
      seriesData.push({
        name : rankRegionName[k],
        data : top[rankRegionName[k]]
      })
    }
    $('#Date_Global_Country_TU').highcharts({
      colors: ["#0000ff", "#ff0066", "#2b908f", "#90ee7e", "#eeaaee",
      "#55BF3B", "#DF5353", "#7798BF", "#aaeeee"],
          chart: {
                  zoomType: 'x',
          backgroundColor: null
          },
          title: {
              text: 'Global Distribution'
          },
          subtitle: {
            text: 'Total Users'
          },
          xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },

              title: {
                  text: 'Date'
              },
              floor: Date.UTC(2015, 00, 01)
          },
          yAxis: {
              title: {
                  text: 'Users'
              },
              min: 0
          },
          tooltip: {
              headerFormat: '<b>{series.name}</b><br>',
              pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
          },

          plotOptions: {
              spline: {
                  marker: {
                      enabled: false
                  }
              }
          },
          credits: {
          enabled: false
        },
          series: seriesData
      });

  },
  Date_TU_China_FX_FC_SC : function(data) {
    var fx = [], fc = [], sc = [];
    for(var i = 0 ; i < data["Date_TU_China_FX"].length; i ++) {
      fx.push([Date.parse(data["Date_TU_China_FX"][i].date), parseInt(data["Date_TU_China_FX"][i]["Total Users"])])
    }

    for(var i = 0 ; i < data["Date_TU_China_FC"].length; i ++) {
      fc.push([Date.parse(data["Date_TU_China_FC"][i].date), parseInt(data["Date_TU_China_FC"][i]["Total Users"])])
    }
    for(var i = 0 ; i < data["Date_TU_China_SC"].length; i ++) {
      sc.push([Date.parse(data["Date_TU_China_SC"][i].date), parseInt(data["Date_TU_China_SC"][i]["Total Users"])])
    }
    $('#Date_TU_China_FX_FC_SC').highcharts({
        colors: ["#0000ff", "#ff0066", "#2b908f", "#90ee7e", "#eeaaee",
      "#55BF3B", "#DF5353", "#7798BF", "#aaeeee"],
          chart: {
              zoomType: 'x',
      backgroundColor: null
          },
          title: {
              text: 'China FX FC SC'
          },
          subtitle: {
              text: 'Total Users'
          },
          xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
              title: {
                text: 'Date'
              },

              floor: Date.UTC(2014, 12, 01)
          },
          yAxis: {
              floor: 0,
              title: {
                  text: 'Users'
              }
          },
          legend: {
              enabled: true
          },
          plotOptions: {
              area: {
                  color: '#0038ff',
                  marker: {
                      enabled: null
                  },
                  lineWidth: 1,
                  states: {
                      hover: {
                          lineWidth: 1
                      }
                  },
                  threshold: null
              }
          },
          credits: {
          enabled: false
        },
          series: [
          {
              name: 'China_FX',
              data: fx
          },
          {
              name: 'China_FC',
              data: fc
          },
          {
              name: 'China_SC',
              data: sc
          }
          ]
      });
  },
  FX_FC_SC_OverlapGrowing : function(data) {
    var fxfc = [], fxsc = [], fcsc = [], fxfcsc = [];
    for(var i = 0 ; i < data["FX_FC_Overlap"].length; i ++) {
      fxfc.push([Date.parse(data["FX_FC_Overlap"][i].date), parseInt(data["FX_FC_Overlap"][i]["Common Users"])])
    }

    for(var i = 0 ; i < data["FX_SC_Overlap"].length; i ++) {
      fxsc.push([Date.parse(data["FX_SC_Overlap"][i].date), parseInt(data["FX_SC_Overlap"][i]["Common Users"])])
    }
    for(var i = 0 ; i < data["FC_SC_Overlap"].length; i ++) {
      fcsc.push([Date.parse(data["FC_SC_Overlap"][i].date), parseInt(data["FC_SC_Overlap"][i]["Common Users"])])
    }
    for(var i = 0 ; i < data["FX_FC_SC_Overlap"].length; i ++) {
      fxfcsc.push([Date.parse(data["FX_FC_SC_Overlap"][i].date), parseInt(data["FX_FC_SC_Overlap"][i]["Common Users"])])
    }
    $('#FX_FC_SC_OverlapGrowing').highcharts({
      colors: ["#0000ff", "#ff0066", "#2b908f", "#90ee7e", "#eeaaee",
    "#55BF3B", "#DF5353", "#7798BF", "#aaeeee"],
        chart: {
      backgroundColor: null,
            zoomType: 'x'
        },
        title: {
            text: 'FX FC SC Users'
        },
        subtitle: {
            text: 'Overlaping Users'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: true
        },
        plotOptions: {
            area: {
                color: '#0038ff',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            name: 'FX&FC',
            data: fxfc
        },{
            name: 'FX&SC',
            data: fxsc
        },
        {
            name: 'FC&SC',
            data: fcsc
        },
        {
            name: 'FX&FC&SC',
            data: fxfcsc
        }]
    });
  },

  VEEN : function(data) {
      (function(venn) {
      "use strict";
      /** given a list of set objects, and their corresponding overlaps.
      updates the (x, y, radius) attribute on each set such that their positions
      roughly correspond to the desired overlaps */
      venn.venn = function(sets, overlaps, parameters) {
          parameters = parameters || {};
          parameters.maxIterations = parameters.maxIterations || 500;
          var lossFunction = parameters.lossFunction || venn.lossFunction;
          var initialLayout = parameters.layoutFunction || venn.greedyLayout;

          // initial layout is done greedily
          sets = initialLayout(sets, overlaps);

          // transform x/y coordinates to a vector to optimize
          var initial = new Array(2*sets.length);
          for (var i = 0; i < sets.length; ++i) {
              initial[2 * i] = sets[i].x;
              initial[2 * i + 1] = sets[i].y;
          }

          // optimize initial layout from our loss function
          var totalFunctionCalls = 0;
          var solution = venn.fmin(
              function(values) {
                  totalFunctionCalls += 1;
                  var current = new Array(sets.length);
                  for (var i = 0; i < sets.length; ++i) {
                      current[i] = {x: values[2 * i],
                                    y: values[2 * i + 1],
                                    radius : sets[i].radius,
                                    size : sets[i].size};
                  }
                  return lossFunction(current, overlaps);
              },
              initial,
              parameters);

          // transform solution vector back to x/y points
          var positions = solution.solution;
          for (i = 0; i < sets.length; ++i) {
              sets[i].x = positions[2 * i];
              sets[i].y = positions[2 * i + 1];
          }

          return sets;
      };

      /** Returns the distance necessary for two circles of radius r1 + r2 to
      have the overlap area 'overlap' */
      venn.distanceFromIntersectArea = function(r1, r2, overlap) {
          // handle complete overlapped circles
          if (Math.min(r1, r2) * Math.min(r1,r2) * Math.PI <= overlap) {
              return Math.abs(r1 - r2);
          }

          return venn.bisect(function(distance) {
              return circleIntersection.circleOverlap(r1, r2, distance) - overlap;
          }, 0, r1 + r2);
      };

      /// gets a matrix of euclidean distances between all sets in venn diagram
      venn.getDistanceMatrix = function(sets, overlaps) {
          // initialize an empty distance matrix between all the points
          var distances = [];
          for (var i = 0; i < sets.length; ++i) {
              distances.push([]);
              for (var j = 0; j < sets.length; ++j) {
                  distances[i].push(0);
              }
          }

          // compute distances between all the points
          for (i = 0; i < overlaps.length; ++i) {
              var current = overlaps[i];
              if (current.sets.length !== 2) {
                  continue;
              }

              var left = current.sets[0],
                  right = current.sets[1],
                  r1 = Math.sqrt(sets[left].size / Math.PI),
                  r2 = Math.sqrt(sets[right].size / Math.PI),
                  distance = venn.distanceFromIntersectArea(r1, r2, current.size);
              distances[left][right] = distances[right][left] = distance;
          }
          return distances;
      };

      /** Lays out a Venn diagram greedily, going from most overlapped sets to
      least overlapped, attempting to position each new set such that the
      overlapping areas to already positioned sets are basically right */
      venn.greedyLayout = function(sets, overlaps) {
          // give each set a default position + radius
          var setOverlaps = {};
          for (var i = 0; i < sets.length; ++i) {
              setOverlaps[i] = [];
              sets[i].radius = Math.sqrt(sets[i].size / Math.PI);
              sets[i].x = sets[i].y = 0;
          }

          // map each set to a list of all the other sets that overlap it
          for (i = 0; i < overlaps.length; ++i) {
              var current = overlaps[i];
              if (current.sets.length !== 2) {
                  continue;
              }

              var left = current.sets[0], right = current.sets[1];
              setOverlaps[left].push ({set:right, size:current.size});
              setOverlaps[right].push({set:left,  size:current.size});
          }

          // get list of most overlapped sets
          var mostOverlapped = [];
          for (var set in setOverlaps) {
              if (setOverlaps.hasOwnProperty(set)) {
                  var size = 0;
                  for (i = 0; i < setOverlaps[set].length; ++i) {
                      size += setOverlaps[set][i].size;
                  }

                  mostOverlapped.push({set: set, size:size});
              }
          }

          // sort by size desc
          function sortOrder(a,b) {
              return b.size - a.size;
          }
          mostOverlapped.sort(sortOrder);

          // keep track of what sets have been laid out
          var positioned = {};
          function isPositioned(element) {
              return element.set in positioned;
          }

          // adds a point to the output
          function positionSet(point, index) {
              sets[index].x = point.x;
              sets[index].y = point.y;
              positioned[index] = true;
          }

          // add most overlapped set at (0,0)
          positionSet({x: 0, y: 0}, mostOverlapped[0].set);

          // get distances between all points
          var distances = venn.getDistanceMatrix(sets, overlaps);

          for (i = 1; i < mostOverlapped.length; ++i) {
              var setIndex = mostOverlapped[i].set,
                  overlap = setOverlaps[setIndex].filter(isPositioned);
              set = sets[setIndex];
              overlap.sort(sortOrder);

              if (overlap.length === 0) {
                  throw "Need overlap information for set " + JSON.stringify( set );
              }

              var points = [];
              for (var j = 0; j < overlap.length; ++j) {
                  // get appropriate distance from most overlapped already added set
                  var p1 = sets[overlap[j].set],
                      d1 = distances[setIndex][overlap[j].set];

                  // sample positions at 90 degrees for maximum aesthetics
                  points.push({x : p1.x + d1, y : p1.y});
                  points.push({x : p1.x - d1, y : p1.y});
                  points.push({y : p1.y + d1, x : p1.x});
                  points.push({y : p1.y - d1, x : p1.x});

                  // if we have at least 2 overlaps, then figure out where the
                  // set should be positioned analytically and try those too
                  for (var k = j + 1; k < overlap.length; ++k) {
                      var p2 = sets[overlap[k].set],
                          d2 = distances[setIndex][overlap[k].set];

                      var extraPoints = circleIntersection.circleCircleIntersection(
                          { x: p1.x, y: p1.y, radius: d1},
                          { x: p2.x, y: p2.y, radius: d2});

                      for (var l = 0; l < extraPoints.length; ++l) {
                          points.push(extraPoints[l]);
                      }
                  }
              }

              // we have some candidate positions for the set, examine loss
              // at each position to figure out where to put it at
              var bestLoss = 1e50, bestPoint = points[0];
              for (j = 0; j < points.length; ++j) {
                  sets[setIndex].x = points[j].x;
                  sets[setIndex].y = points[j].y;
                  var loss = venn.lossFunction(sets, overlaps);
                  if (loss < bestLoss) {
                      bestLoss = loss;
                      bestPoint = points[j];
                  }
              }

              positionSet(bestPoint, setIndex);
          }

          return sets;
      };

      /// Uses multidimensional scaling to approximate a first layout here
      venn.classicMDSLayout = function(sets, overlaps) {
          // get the distance matrix
          var distances = venn.getDistanceMatrix(sets, overlaps);

          // get positions for each set
          var positions = mds.classic(distances);

          // translate back to (x,y,radius) coordinates
          for (var i = 0; i < sets.length; ++i) {
              sets[i].x = positions[i][0];
              sets[i].y = positions[i][1];
              sets[i].radius = Math.sqrt(sets[i].size / Math.PI);
          }
          return sets;
      };

      /** Given a bunch of sets, and the desired overlaps between these sets - computes
      the distance from the actual overlaps to the desired overlaps. Note that
      this method ignores overlaps of more than 2 circles */
      venn.lossFunction = function(sets, overlaps) {
          var output = 0;

          function getCircles(indices) {
              return indices.map(function(i) { return sets[i]; });
          }

          for (var i = 0; i < overlaps.length; ++i) {
              var area = overlaps[i], overlap;
              if (area.sets.length == 2) {
                  var left = sets[area.sets[0]],
                      right = sets[area.sets[1]];
                  overlap = circleIntersection.circleOverlap(left.radius, right.radius,
                                  circleIntersection.distance(left, right));
              } else {
                  overlap = circleIntersection.intersectionArea(getCircles(area.sets));
              }

              output += (overlap - area.size) * (overlap - area.size);
          }

          return output;
      };

      /** Scales a solution from venn.venn or venn.greedyLayout such that it fits in
      a rectangle of width/height - with padding around the borders. */
      venn.scaleSolution = function(solution, width, height, padding) {
          var minMax = function(d) {
              var hi = Math.max.apply(null, solution.map(
                                      function(c) { return c[d] + c.radius; } )),
                  lo = Math.min.apply(null, solution.map(
                                      function(c) { return c[d] - c.radius;} ));
              return {max:hi, min:lo};
          };

          width -= 2*padding;
          height -= 2*padding;

          var xRange = minMax('x'),
              yRange = minMax('y'),
              xScaling = width  / (xRange.max - xRange.min),
              yScaling = height / (yRange.max - yRange.min),
              scaling = Math.min(yScaling, xScaling);

          for (var i = 0; i < solution.length; ++i) {
              var set = solution[i];
              set.radius = scaling * set.radius;
              set.x = padding + (set.x - xRange.min) * scaling;
              set.y = padding + (set.y - yRange.min) * scaling;
          }
          solution.scaling = scaling;

          return solution;
      };

      function weightedSum(a, b) {
          var ret = new Array(a[1].length || 0);
          for (var j = 0; j < ret.length; ++j) {
              ret[j] = a[0] * a[1][j] + b[0] * b[1][j];
          }
          return ret;
      }

      function centerVennDiagram( diagram, width, height, padding ) {
          var diagramBoundaries;
          var allowedWidth = width - ( 2 * ( padding || 0 ) );
          var allowedHeight = height - ( 2 * ( padding || 0 ) );
          var scale;
          var transformX, transformY;
          var transform = "";

          if ( diagram ) {
              diagramBoundaries = diagram[ 0 ][ 0 ].getBBox();
              if ( diagramBoundaries && width && height ) {

                  //  See if we need to scale to fit the width/height
                  if ( diagramBoundaries.width > allowedWidth ) {
                      scale = allowedWidth / diagramBoundaries.width;
                  }
                  if ( diagramBoundaries.height > allowedHeight ) {
                      if ( !scale || ( allowedHeight / diagramBoundaries.height ) < scale ) {
                          scale = allowedHeight / diagramBoundaries.height;
                      }
                  }

                  if ( scale ) {
                      transform = "scale(" + scale + ")";
                  }
                  else {
                      scale = 1;
                  }

                  transformX = Math.floor( ( allowedWidth - ( diagramBoundaries.width * scale ) ) / 2 );
                  transformY = Math.floor( ( allowedHeight - ( diagramBoundaries.height * scale ) ) / 2 );
                  diagram.attr( "transform", "translate(" + transformX + ","  + transformY + ") " + transform );
              }
          }
      }

      /** finds the zeros of a function, given two starting points (which must
       * have opposite signs */
      venn.bisect = function(f, a, b, parameters) {
          parameters = parameters || {};
          var maxIterations = parameters.maxIterations || 100,
              tolerance = parameters.tolerance || 1e-10,
              fA = f(a),
              fB = f(b),
              delta = b - a;

          if (fA * fB > 0) {
              throw "Initial bisect points must have opposite signs";
          }

          if (fA === 0) return a;
          if (fB === 0) return b;

          for (var i = 0; i < maxIterations; ++i) {
              delta /= 2;
              var mid = a + delta,
                  fMid = f(mid);

              if (fMid * fA >= 0) {
                  a = mid;
              }

              if ((Math.abs(delta) < tolerance) || (fMid === 0)) {
                  return mid;
              }
          }
          return a + delta;
      };

      /** minimizes a function using the downhill simplex method */
      venn.fmin = function(f, x0, parameters) {
          parameters = parameters || {};

          var maxIterations = parameters.maxIterations || x0.length * 200,
              nonZeroDelta = parameters.nonZeroDelta || 1.1,
              zeroDelta = parameters.zeroDelta || 0.001,
              minErrorDelta = parameters.minErrorDelta || 1e-5,
              rho = parameters.rho || 1,
              chi = parameters.chi || 2,
              psi = parameters.psi || -0.5,
              sigma = parameters.sigma || 0.5,
              callback = parameters.callback;

          // initialize simplex.
          var N = x0.length,
              simplex = new Array(N + 1);
          simplex[0] = x0;
          simplex[0].fx = f(x0);
          for (var i = 0; i < N; ++i) {
              var point = x0.slice();
              point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
              simplex[i+1] = point;
              simplex[i+1].fx = f(point);
          }

          var sortOrder = function(a, b) { return a.fx - b.fx; };

          for (var iteration = 0; iteration < maxIterations; ++iteration) {
              simplex.sort(sortOrder);
              if (callback) {
                  callback(simplex);
              }

              if (Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) {
                  break;
              }

              // compute the centroid of all but the worst point in the simplex
              var centroid = new Array(N);
              for (i = 0; i < N; ++i) {
                  centroid[i] = 0;
                  for (var j = 0; j < N; ++j) {
                      centroid[i] += simplex[j][i];
                  }
                  centroid[i] /= N;
              }

              // reflect the worst point past the centroid  and compute loss at reflected
              // point
              var worst = simplex[N];
              var reflected = weightedSum([1+rho, centroid], [-rho, worst]);
              reflected.fx = f(reflected);

              var replacement = reflected;

              // if the reflected point is the best seen, then possibly expand
              if (reflected.fx <= simplex[0].fx) {
                  var expanded = weightedSum([1+chi, centroid], [-chi, worst]);
                  expanded.fx = f(expanded);
                  if (expanded.fx < reflected.fx) {
                      replacement = expanded;
                  }
              }

              // if the reflected point is worse than the second worst, we need to
              // contract
              else if (reflected.fx >= simplex[N-1].fx) {
                  var shouldReduce = false;
                  var contracted;

                  if (reflected.fx <= worst.fx) {
                      // do an inside contraction
                      contracted = weightedSum([1+psi, centroid], [-psi, worst]);
                      contracted.fx = f(contracted);
                      if (contracted.fx < worst.fx) {
                          replacement = contracted;
                      } else {
                          shouldReduce = true;
                      }
                  } else {
                      // do an outside contraction
                      contracted = weightedSum([1-psi * rho, centroid], [psi*rho, worst]);
                      contracted.fx = f(contracted);
                      if (contracted.fx <= reflected.fx) {
                          replacement = contracted;
                      } else {
                          shouldReduce = true;
                      }
                  }

                  if (shouldReduce) {
                      // do reduction. doesn't actually happen that often
                      for (i = 1; i < simplex.length; ++i) {
                          simplex[i] = weightedSum([1 - sigma, simplex[0]],
                                                   [sigma - 1, simplex[i]]);
                          simplex[i].fx = f(simplex[i]);
                      }
                  }
              }

              simplex[N] = replacement;
          }

          simplex.sort(sortOrder);
          return {f : simplex[0].fx,
                  solution : simplex[0]};
      };

      /** returns a svg path of the intersection area of a bunch of circles */
      venn.intersectionAreaPath = function(circles) {
          var stats = {};
          circleIntersection.intersectionArea(circles, stats);
          var arcs = stats.arcs;

          if (arcs.length == 0) {
              return "M 0 0";
          }

          var ret = ["\nM", arcs[0].p2.x, arcs[0].p2.y];
          for (var i = 0; i < arcs.length; ++i) {
              var arc = arcs[i], r = arc.circle.radius, wide = arc.width > r;
              ret.push("\nA", r, r, 0, wide ? 1 : 0, 1, arc.p1.x, arc.p1.y);
          }

          return ret.join(" ");
      }

      venn.drawD3Diagram = function(element, dataset, width, height, parameters) {
          parameters = parameters || {};

          var colours = d3.scale.category10(),
              circleFillColours = parameters.circleFillColours || colours,
              circleStrokeColours = parameters.circleStrokeColours || circleFillColours,
              circleStrokeWidth = parameters.circleStrokeWidth || function(i) { return 0; },
              textFillColours = parameters.textFillColours || colours,
              textStrokeColours = parameters.textStrokeColours || textFillColours,
              nodeOpacity = parameters.opacity || 0.3,
              padding = parameters.padding || 6;

          dataset = venn.scaleSolution(dataset, width, height, padding);
          var svg = element.append("svg")
                  .attr("width", width)
                  .attr("height", height);

          var diagram = svg.append( "g" );

          var nodes = diagram.append("g").selectAll("circle")
                           .data(dataset)
                           .enter()
                           .append("g");

          var circles = nodes.append("circle")
                 .attr("r",  function(d) { return d.radius; })
                 .style("fill-opacity", nodeOpacity)
                 .attr("cx", function(d) { return d.x; })
                 .attr("cy", function(d) { return d.y; })
                 .style("stroke", function(d, i) { return circleStrokeColours(i); })
                 .style("stroke-width", function(d, i) { return circleStrokeWidth(i); })
                 .style("fill", function(d, i) { return circleFillColours(i); });

          var text = nodes.append("text")
                 .attr("x", function(d) { return d.x; })
                 .attr("y", function(d) { return d.y; })
                 .attr("text-anchor", "middle")
                 .attr("dy", "0.35em")
                 .style("stroke", function(d, i) { return textStrokeColours(i); })
                 .style("fill", function(d, i) { return textFillColours(i); })
                 .text(function(d) { return d.label; });

          centerVennDiagram( diagram, width, height, padding );

          return {'svg' : svg,
                  'nodes' : nodes,
                  'circles' : circles,
                  'text' : text };
      };

      venn.updateD3Diagram = function(element, dataset) {
          var svg = element.select("svg"),
              width = parseInt(svg.attr('width'), 10),
              height = parseInt(svg.attr('height'), 10);

          dataset = venn.scaleSolution(dataset, width, height, 6);
          element.selectAll("circle")
              .data(dataset)
              .transition()
              .duration(400)
              .attr("cx", function(d) { return d.x; })
              .attr("cy", function(d) { return d.y; })
              .attr("r",  function(d) { return d.radius; });

          element.selectAll("text")
              .data(dataset)
              .transition()
              .duration(400)
              .text(function(d) { return d.label; })
              .attr("x", function(d) { return d.x; })
              .attr("y", function(d) { return d.y; });
      };
  }(window.venn = window.venn || {}));
  (function(circleIntersection) {
      "use strict";
      var SMALL = 1e-10;

      /** Returns the intersection area of a bunch of circles (where each circle
       is an object having an x,y and radius property) */
      circleIntersection.intersectionArea = function(circles, stats) {
          // get all the intersection points of the circles
          var intersectionPoints = getIntersectionPoints(circles);

          // filter out points that aren't included in all the circles
          var innerPoints = intersectionPoints.filter(function (p) {
              return circleIntersection.containedInCircles(p, circles);
          });

          var arcArea = 0, polygonArea = 0, arcs = [], i;

          // if we have intersection points that are within all the circles,
          // then figure out the area contained by them
          if (innerPoints.length > 1) {
              // sort the points by angle from the center of the polygon, which lets
              // us just iterate over points to get the edges
              var center = circleIntersection.getCenter(innerPoints);
              for (i = 0; i < innerPoints.length; ++i ) {
                  var p = innerPoints[i];
                  p.angle = Math.atan2(p.x - center.x, p.y - center.y);
              }
              innerPoints.sort(function(a,b) { return b.angle - a.angle;});

              // iterate over all points, get arc between the points
              // and update the areas
              var p2 = innerPoints[innerPoints.length - 1];
              for (i = 0; i < innerPoints.length; ++i) {
                  var p1 = innerPoints[i];

                  // polygon area updates easily ...
                  polygonArea += (p2.x + p1.x) * (p1.y - p2.y);

                  // updating the arc area is a little more involved
                  var midPoint = {x : (p1.x + p2.x) / 2,
                                  y : (p1.y + p2.y) / 2},
                      arc = null;

                  for (var j = 0; j < p1.parentIndex.length; ++j) {
                      if (p2.parentIndex.indexOf(p1.parentIndex[j]) > -1) {
                          // figure out the angle halfway between the two points
                          // on the current circle
                          var circle = circles[p1.parentIndex[j]],
                              a1 = Math.atan2(p1.x - circle.x, p1.y - circle.y),
                              a2 = Math.atan2(p2.x - circle.x, p2.y - circle.y);

                          var angleDiff = (a2 - a1);
                          if (angleDiff < 0) {
                              angleDiff += 2*Math.PI;
                          }

                          // and use that angle to figure out the width of the
                          // arc
                          var a = a2 - angleDiff/2,
                              width = circleIntersection.distance(midPoint, {
                                  x : circle.x + circle.radius * Math.sin(a),
                                  y : circle.y + circle.radius * Math.cos(a)
                              });

                          // pick the circle whose arc has the smallest width
                          if ((arc === null) || (arc.width > width)) {
                              arc = { circle : circle,
                                      width : width,
                                      p1 : p1,
                                      p2 : p2};
                          }
                      }
                  }
                  arcs.push(arc);
                  arcArea += circleIntersection.circleArea(arc.circle.radius, arc.width);
                  p2 = p1;
              }
          } else {
              // no intersection points, is either disjoint - or is completely
              // overlapped. figure out which by examining the smallest circle
              var smallest = circles[0];
              for (i = 1; i < circles.length; ++i) {
                  if (circles[i].radius < smallest.radius) {
                      smallest = circles[i];
                  }
              }

              // make sure the smallest circle is completely contained in all
              // the other circles
              var disjoint = false;
              for (i = 0; i < circles.length; ++i) {
                  if (circleIntersection.distance(circles[i], smallest) > Math.abs(smallest.radius - circles[i].radius)) {
                      disjoint = true;
                      break;
                  }
              }

              if (disjoint) {
                  arcArea = polygonArea = 0;

              } else {
                  arcArea = smallest.radius * smallest.radius * Math.PI;
                  arcs.push({circle : smallest,
                             p1: { x: smallest.x,        y : smallest.y + smallest.radius},
                             p2: { x: smallest.x - SMALL, y : smallest.y + smallest.radius},
                             width : smallest.radius * 2 });
              }
          }

          polygonArea /= 2;
          if (stats) {
              stats.area = arcArea + polygonArea;
              stats.arcArea = arcArea;
              stats.polygonArea = polygonArea;
              stats.arcs = arcs;
              stats.innerPoints = innerPoints;
              stats.intersectionPoints = intersectionPoints;
          }

          return arcArea + polygonArea;
      };

      /** returns whether a point is contained by all of a list of circles */
      circleIntersection.containedInCircles = function(point, circles) {
          for (var i = 0; i < circles.length; ++i) {
              if (circleIntersection.distance(point, circles[i]) > circles[i].radius + SMALL) {
                  return false;
              }
          }
          return true;
      };

      /** Gets all intersection points between a bunch of circles */
      function getIntersectionPoints(circles) {
          var ret = [];
          for (var i = 0; i < circles.length; ++i) {
              for (var j = i + 1; j < circles.length; ++j) {
                  var intersect = circleIntersection.circleCircleIntersection(circles[i],
                                                                circles[j]);
                  for (var k = 0; k < intersect.length; ++k) {
                      var p = intersect[k];
                      p.parentIndex = [i,j];
                      ret.push(p);
                  }
              }
          }
          return ret;
      }

      circleIntersection.circleIntegral = function(r, x) {
          var y = Math.sqrt(r * r - x * x);
          return x * y + r * r * Math.atan2(x, y);
      };

      /** Returns the area of a circle of radius r - up to width */
      circleIntersection.circleArea = function(r, width) {
          return circleIntersection.circleIntegral(r, width - r) - circleIntersection.circleIntegral(r, -r);
      };


      /** euclidean distance between two points */
      circleIntersection.distance = function(p1, p2) {
          return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                           (p1.y - p2.y) * (p1.y - p2.y));
      };


      /** Returns the overlap area of two circles of radius r1 and r2 - that
      have their centers separated by distance d. Simpler faster
      circle intersection for only two circles */
      circleIntersection.circleOverlap = function(r1, r2, d) {
          // no overlap
          if (d >= r1 + r2) {
              return 0;
          }

          // completely overlapped
          if (d <= Math.abs(r1 - r2)) {
              return Math.PI * Math.min(r1, r2) * Math.min(r1, r2);
          }

          var w1 = r1 - (d * d - r2 * r2 + r1 * r1) / (2 * d),
              w2 = r2 - (d * d - r1 * r1 + r2 * r2) / (2 * d);
          return circleIntersection.circleArea(r1, w1) + circleIntersection.circleArea(r2, w2);
      };


      /** Given two circles (containing a x/y/radius attributes),
      returns the intersecting points if possible.
      note: doesn't handle cases where there are infinitely many
      intersection points (circles are equivalent):, or only one intersection point*/
      circleIntersection.circleCircleIntersection = function(p1, p2) {
          var d = circleIntersection.distance(p1, p2),
              r1 = p1.radius,
              r2 = p2.radius;

          // if to far away, or self contained - can't be done
          if ((d >= (r1 + r2)) || (d <= Math.abs(r1 - r2))) {
              return [];
          }

          var a = (r1 * r1 - r2 * r2 + d * d) / (2 * d),
              h = Math.sqrt(r1 * r1 - a * a),
              x0 = p1.x + a * (p2.x - p1.x) / d,
              y0 = p1.y + a * (p2.y - p1.y) / d,
              rx = -(p2.y - p1.y) * (h / d),
              ry = -(p2.x - p1.x) * (h / d);

          return [{ x: x0 + rx, y : y0 - ry },
                  { x: x0 - rx, y : y0 + ry }];
      };

      /** Returns the center of a bunch of points */
      circleIntersection.getCenter = function(points) {
          var center = { x: 0, y: 0};
          for (var i =0; i < points.length; ++i ) {
              center.x += points[i].x;
              center.y += points[i].y;
          }
          center.x /= points.length;
          center.y /= points.length;
          return center;
      };
  }(window.circleIntersection = window.circleIntersection || {}));

  var sets = [
         {"label": 'China_FX', "size": data["China_FX_FC_SC_TU"][0]["China_FX"]},
         {"label": "China_FC", "size": data["China_FX_FC_SC_TU"][0]["China_FC"]},
         {"label": "China_SC", "size": data["China_FX_FC_SC_TU"][0]["China_SC"]}
    ],

    overlaps = [
         {"sets": [0, 1], "size": data["China_FX_FC_SC_TU"][0]["China_FX_FC"]},
         {"sets": [0, 2], "size": data["China_FX_FC_SC_TU"][0]["China_FX_SC"]},
         {"sets": [1, 2], "size": data["China_FX_FC_SC_TU"][0]["China_FC_SC"]},
         {"sets": [0, 1, 2], "size": data["China_FX_FC_SC_TU"][0]["China_FX_FC_SC"] },
       ],
  // get positions for each set
  sets = venn.venn(sets, overlaps);

  // draw the diagram in the 'venn' div
  var diagram = venn.drawD3Diagram(d3.select(".venn"), sets, 500, 500);

  // add a tooltip showing the size of each set/intersection
  var tooltip = d3.select("body").append("div")
      .attr("class", "venntooltip");

  d3.selection.prototype.moveParentToFront = function() {
    return this.each(function(){
      this.parentNode.parentNode.appendChild(this.parentNode);
    });
  };

  // hover on all the circles
  diagram.circles
      .style("stroke-opacity", 0)
      .style("stroke", "white")
      .style("stroke-width", "2")
      .on("mousemove", function() {
          tooltip.style("left", (d3.event.pageX) + "px")
                 .style("top", (d3.event.pageY - 28) + "px");
      })
      .on("mouseover", function(d, i) {
          var selection = d3.select(this);
          d3.select(this).moveParentToFront()
              .transition()
              .style("fill-opacity", .5)
              .style("stroke-opacity", 1);

          tooltip.transition().style("opacity", .9);
          tooltip.text(d.size + " users");
      })
      .on("mouseout", function(d, i) {
          d3.select(this).transition()
              .style("fill-opacity", .3)
              .style("stroke-opacity", 0);
          tooltip.transition().style("opacity", 0);
      });

  // draw a path around each intersection area, add hover there as well
  diagram.svg.select("g").selectAll("path")
      .data(overlaps)
      .enter()
      .append("path")
      .attr("d", function(d) {
          return venn.intersectionAreaPath(d.sets.map(function(j) { return sets[j]; }));
      })
      .style("fill-opacity","0")
      .style("fill", "black")
      .style("stroke-opacity", 0)
      .style("stroke", "white")
      .style("stroke-width", "2")
      .on("mouseover", function(d, i) {
          d3.select(this).transition()
              .style("fill-opacity", .1)
              .style("stroke-opacity", 1);
          tooltip.transition().style("opacity", .9);
          tooltip.text(d.size + " users");
      })
      .on("mouseout", function(d, i) {
          d3.select(this).transition()
              .style("fill-opacity", 0)
              .style("stroke-opacity", 0);
          tooltip.transition().style("opacity", 0);
      })
      .on("mousemove", function() {
          tooltip.style("left", (d3.event.pageX) + "px")
                 .style("top", (d3.event.pageY - 28) + "px");
      })

  },
  Date_NU_Global : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Date_NU_Global"].length; i ++) {
      dng.push([Date.parse(data["Date_NU_Global"][i].date), parseInt(data["Date_NU_Global"][i]["New Users"])])
    }
    $('#Date_NU_Global').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Daily New Users'
        },
        subtitle: {
            text: 'Global'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#07a86d',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'New Users',
            data: dng
        }]
    });
  },
  Date_NU_China : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Date_NU_China"].length; i ++) {
      dng.push([Date.parse(data["Date_NU_China"][i].date), parseInt(data["Date_NU_China"][i]["New Users"])])
    }
    $('#Date_NU_China').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Daily New Users'
        },
        subtitle: {
            text: 'China'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#07a86d',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'New Users',
            data: dng
        }]
    });
  },
  Week_NU_Global : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Week_NU_Global"].length; i ++) {
      dng.push([Date.parse(data["Week_NU_Global"][i].date), parseInt(data["Week_NU_Global"][i]["New Users"])])
    }
    $('#Week_NU_Global').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Weekly New Users'
        },
        subtitle: {
            text: 'Global'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#07a86d',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'New Users',
            data: dng
        }]
    });
  },
  Week_NU_China : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Week_NU_China"].length; i ++) {
      dng.push([Date.parse(data["Week_NU_China"][i].date), parseInt(data["Week_NU_China"][i]["New Users"])])
    }
    $('#Week_NU_China').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Weekly New Users'
        },
        subtitle: {
            text: 'China'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#07a86d',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'New Users',
            data: dng
        }]
    });
  },
  Month_NU_Global : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Month_NU_Global"].length; i ++) {
      dng.push([Date.parse(data["Month_NU_Global"][i].date), parseInt(data["Month_NU_Global"][i]["New Users"])])
    }
    $('#Month_NU_Global').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Monthly New Users'
        },
        subtitle: {
            text: 'Global'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#07a86d',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'New Users',
            data: dng
        }]
    });
  },
  Month_NU_China : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Month_NU_China"].length; i ++) {
      dng.push([Date.parse(data["Month_NU_China"][i].date), parseInt(data["Month_NU_China"][i]["New Users"])])
    }
    $('#Month_NU_China').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Monthly New Users'
        },
        subtitle: {
            text: 'China'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#07a86d',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'New Users',
            data: dng
        }]
    });
  },

Date_AU_Global_FX : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Date_AU_Global_FX"].length; i ++) {
      dng.push([Date.parse(data["Date_AU_Global_FX"][i].date), parseInt(data["Date_AU_Global_FX"][i]["Active Users"])])
    }
    $('#Date_AU_Global_FX').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
            text: 'Global_FX'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2014, 12, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },

        credits: {
          enabled: false
        },

        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Date_AU_China_FX : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Date_AU_China_FX"].length; i ++) {
      dng.push([Date.parse(data["Date_AU_China_FX"][i].date), parseInt(data["Date_AU_China_FX"][i]["Active Users"])])
    }
    $('#Date_AU_China_FX').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
            text: 'China_FX'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2014, 12, 01)
        },
        yAxis: {
            ceiling: 8000,
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },

        credits: {
          enabled: false
        },

        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Date_AU_China_FC : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Date_AU_China_FC"].length; i ++) {
      dng.push([Date.parse(data["Date_AU_China_FC"][i].date), parseInt(data["Date_AU_China_FC"][i]["Active Users"])])
    }
    $('#Date_AU_China_FC').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
            text: 'China_FC'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 04, 03)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },

        credits: {
          enabled: false
        },

        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },

Date_AU_China_SC : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Date_AU_China_SC"].length; i ++) {
      dng.push([Date.parse(data["Date_AU_China_SC"][i].date), parseInt(data["Date_AU_China_SC"][i]["Active Users"])])
    }
    $('#Date_AU_China_SC').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
            text: 'China_SC'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 06, 06)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },

        credits: {
          enabled: false
        },

        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Week_AU_Global_FX : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Week_AU_Global_FX"].length; i ++) {
      dng.push([Date.parse(data["Week_AU_Global_FX"][i].date), parseInt(data["Week_AU_Global_FX"][i]["Active Users"])])
    }
    $('#Week_AU_Global_FX').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
            text: 'Global_FX'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2014, 12, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Week_AU_China_FX : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Week_AU_China_FX"].length; i ++) {
      dng.push([Date.parse(data["Week_AU_China_FX"][i].date), parseInt(data["Week_AU_China_FX"][i]["Active Users"])])
    }
    $('#Week_AU_China_FX').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
            text: 'China_FX'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2014, 12, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Week_AU_China_FC : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Week_AU_China_FC"].length; i ++) {
      dng.push([Date.parse(data["Week_AU_China_FC"][i].date), parseInt(data["Week_AU_China_FC"][i]["Active Users"])])
    }
    $('#Week_AU_China_FC').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
            text: 'China_FC'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 04, 04)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Week_AU_China_SC : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Week_AU_China_SC"].length; i ++) {
      dng.push([Date.parse(data["Week_AU_China_SC"][i].date), parseInt(data["Week_AU_China_SC"][i]["Active Users"])])
    }
    $('#Week_AU_China_SC').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
            text: 'China_SC'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 06, 06)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Month_AU_Global_FX : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Month_AU_Global_FX"].length; i ++) {
      dng.push([Date.parse(data["Month_AU_Global_FX"][i].date), parseInt(data["Month_AU_Global_FX"][i]["Active Users"])])
    }
    $('#Month_AU_Global_FX').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
            text: 'Global_FX'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2014, 12, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Month_AU_China_FX : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Month_AU_China_FX"].length; i ++) {
      dng.push([Date.parse(data["Month_AU_China_FX"][i].date), parseInt(data["Month_AU_China_FX"][i]["Active Users"])])
    }
    $('#Month_AU_China_FX').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
            text: 'China_FX'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2014, 12, 01)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Month_AU_China_FC : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Month_AU_China_FC"].length; i ++) {
      dng.push([Date.parse(data["Month_AU_China_FC"][i].date), parseInt(data["Month_AU_China_FC"][i]["Active Users"])])
    }
    $('#Month_AU_China_FC').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
            text: 'China_FC'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 04, 04)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
Month_AU_China_SC : function(data) {
    var dng = [];
    for(var i = 0 ; i < data["Month_AU_China_SC"].length; i ++) {
      dng.push([Date.parse(data["Month_AU_China_SC"][i].date), parseInt(data["Month_AU_China_SC"][i]["Active Users"])])
    }
    $('#Month_AU_China_SC').highcharts({
        chart: {
            type: 'spline',
            zoomType: 'x',
      backgroundColor: null
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
            text: 'China_SC'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
              text: 'Date'
            },
            floor: Date.UTC(2015, 06, 06)
        },
        yAxis: {
            floor: 0,
            title: {
                text: 'Users'
            }
        },
        legend: {
            enabled: false
        },
        plotOptions: {
            area: {
                color: '#f9105e',
                marker: {
                    enabled: null
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null
            }
        },
        credits: {
          enabled: false
        },
        series: [{
            type: 'area',
            name: 'Active Users',
            data: dng
        }]
    });
  },
  Date_School_TU_China: function(data) {
    var dstc = [];
    for(var k in data["Top10_School_China"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_School_TU_China"].length; i ++) {
        if(data["Date_School_TU_China"][i].school == data["Top10_School_China"][0][k])
          dataArr.push([Date.parse(data["Date_School_TU_China"][i].date), parseInt(data["Date_School_TU_China"][i]["Total Users"])])
      }
      dstc.push({
        name : data["Top10_School_China"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_TU_China').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'China School Top10'
        },
        subtitle: {
          text: 'Total Users'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dstc
    });
  },

  Date_School_TU_Global: function(data) {
    var dstg = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_School_TU_Global"].length; i ++) {
        if(data["Date_School_TU_Global"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([Date.parse(data["Date_School_TU_Global"][i].date), parseInt(data["Date_School_TU_Global"][i]["Total Users"])])
      }
      dstg.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_TU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Global School Top10'
        },
        subtitle: {
          text: 'Total Users'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dstg
    });
  },

  Date_School_NU_Global: function(data) {
    var dsng = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_School_NU_Global"].length; i ++) {
        if(data["Date_School_NU_Global"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([Date.parse(data["Date_School_NU_Global"][i].date), parseInt(data["Date_School_NU_Global"][i]["New Users"])])
      }
      dsng.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_NU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily New Users'
        },
        subtitle: {
          text: 'Global Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dsng
    });
  },

  Date_School_NU_China: function(data) {
    var dsng = [];
    for(var k in data["Top10_School_China"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_School_NU_China"].length; i ++) {
        if(data["Date_School_NU_China"][i].school == data["Top10_School_China"][0][k])
          dataArr.push([Date.parse(data["Date_School_NU_China"][i].date), parseInt(data["Date_School_NU_China"][i]["New Users"])])
      }
      dsng.push({
        name : data["Top10_School_China"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_NU_China').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily New Users'
        },
        subtitle: {
          text: 'China Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dsng
    });
  },

  Week_School_NU_Global: function(data) {
    var wsng = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = []
      for(var i = 0; i < data["Week_School_NU_Global"].length; i++) {
        if(data["Week_School_NU_Global"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([
            Date.parse(data["Week_School_NU_Global"][i].date),
            parseInt(data["Week_School_NU_Global"][i]["New Users"])
            ])
      }
      wsng.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }

    $('#Week_School_NU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly New Users'
        },
        subtitle: {
          text: 'Global Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsng
    });
  },

  Week_School_NU_China: function(data) {
    var wsnc = [];
    for(var k in data["Top10_School_China"][0]) {
      var dataArr = []
      for(var i = 0; i < data["Week_School_NU_China"].length; i++) {
        if(data["Week_School_NU_China"][i].school == data["Top10_School_China"][0][k])
          dataArr.push([
            Date.parse(data["Week_School_NU_China"][i].date),
            parseInt(data["Week_School_NU_China"][i]["New Users"])
            ])
      }
      wsnc.push({
        name : data["Top10_School_China"][0][k],
        data : dataArr
      })
      // inner loop: for every item in Week_School_NU_China, see if it is in the 1st of the Top10_School_China,
      // and if it is, dataArr.push([ ... ])
      // outer loop: repeat this for all 10 items in Top10_School_China, and it if is, wsnc.push({ ... })
    }
    $('#Week_School_NU_China').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly New Users'
        },
        subtitle: {
          text: 'China Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsnc
    });
  },

  Month_School_NU_Global: function(data) {
    var msng = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = []
      for(var i = 0; i < data["Month_School_NU_Global"].length; i++) {
        if(data["Month_School_NU_Global"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([
            Date.parse(data["Month_School_NU_Global"][i].date),
            parseInt(data["Month_School_NU_Global"][i]["New Users"])
            ])
      }
      msng.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_NU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly New Users'
        },
        subtitle: {
          text: 'Global Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msng
    });
  },

  Month_School_NU_China: function(data) {
    // inner loop: for every item in Week_School_NU_China, see if it is in the 1st of the Top10_School_China,
    // and if it is, dataArr.push([ ... ])
    // outer loop: repeat this for all 10 items in Top10_School_China, and it if is, wsnc.push({ ... })
    var msnc = [];
    for(var k in data["Top10_School_China"][0]) {
      var dataArr = []
      for(var i = 0; i < data["Month_School_NU_China"].length; i++) {
        if(data["Month_School_NU_China"][i].school == data["Top10_School_China"][0][k])
          dataArr.push([
            Date.parse(data["Month_School_NU_China"][i].date),
            parseInt(data["Month_School_NU_China"][i]["New Users"])
            ])
      }
      msnc.push({
        name : data["Top10_School_China"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_NU_China').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly New Users'
        },
        subtitle: {
          text: 'China Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msnc
    });
  },

  //Date_School_AU_Global_FX
  Date_School_AU_Global_FX: function(data) {
    var dsag = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Date_School_AU_Global_FX"].length; i ++) {
        if(data["Date_School_AU_Global_FX"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([Date.parse(data["Date_School_AU_Global_FX"][i].date), parseInt(data["Date_School_AU_Global_FX"][i]["Active Users"])])
      }
      dsag.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_AU_Global_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
          text: 'Global FX Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dsag
    });
  },
  // Week_School_AU_Global_FX
  Week_School_AU_Global_FX: function(data) {
    var wsag = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Week_School_AU_Global_FX"].length; i ++) {
        if(data["Week_School_AU_Global_FX"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([Date.parse(data["Week_School_AU_Global_FX"][i].date), parseInt(data["Week_School_AU_Global_FX"][i]["Active Users"])])
      }
      wsag.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }
    $('#Week_School_AU_Global_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
          text: 'Global FX Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsag
    });
  },

  // Month_School_AU_Global_FX
  Month_School_AU_Global_FX: function(data) {
    var msag = [];
    for(var k in data["Top10_School_Global"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Month_School_AU_Global_FX"].length; i ++) {
        if(data["Month_School_AU_Global_FX"][i].school == data["Top10_School_Global"][0][k])
          dataArr.push([Date.parse(data["Month_School_AU_Global_FX"][i].date), parseInt(data["Month_School_AU_Global_FX"][i]["Active Users"])])
      }
      msag.push({
        name : data["Top10_School_Global"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_AU_Global_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
          text: 'Global FX Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msag
    });
  },

  //Date_School_AU_China_FX
  Date_School_AU_China_FX: function(data) {
    var dsacfx = [];
    for(var k in data["Top10_School_China_FX"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Date_School_AU_China_FX"].length; i ++) {
        if(data["Date_School_AU_China_FX"][i].school == data["Top10_School_China_FX"][0][k])
          dataArr.push([Date.parse(data["Date_School_AU_China_FX"][i].date), parseInt(data["Date_School_AU_China_FX"][i]["Active Users"])])
      }
      dsacfx.push({
        name : data["Top10_School_China_FX"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_AU_China_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
          text: 'China FX Top10 School'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dsacfx
    });
  },
  // Week_School_AU_China_FX
  Week_School_AU_China_FX: function(data) {
    var wsacfx = [];
    for(var k in data["Top10_School_China_FX"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Week_School_AU_China_FX"].length; i ++) {
        if(data["Week_School_AU_China_FX"][i].school == data["Top10_School_China_FX"][0][k])
          dataArr.push([Date.parse(data["Week_School_AU_China_FX"][i].date), parseInt(data["Week_School_AU_China_FX"][i]["Active Users"])])
      }
      wsacfx.push({
        name : data["Top10_School_China_FX"][0][k],
        data : dataArr
      })
    }
    $('#Week_School_AU_China_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
          text: 'China FX Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 05, 20)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsacfx
    });
  },

  // Month_School_AU_China_FX
  Month_School_AU_China_FX: function(data) {
    var msacfx = [];
    for(var k in data["Top10_School_China_FX"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Month_School_AU_China_FX"].length; i ++) {
        if(data["Month_School_AU_China_FX"][i].school == data["Top10_School_China_FX"][0][k])
          dataArr.push([Date.parse(data["Month_School_AU_China_FX"][i].date), parseInt(data["Month_School_AU_China_FX"][i]["Active Users"])])
      }
      msacfx.push({
        name : data["Top10_School_China_FX"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_AU_China_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
          text: 'China FX Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 06, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msacfx
    });
  },

  //Date_School_AU_China_FC
  Date_School_AU_China_FC: function(data) {
    var dsacfc = [];
    for(var k in data["Top10_School_China_FC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Date_School_AU_China_FC"].length; i ++) {
        if(data["Date_School_AU_China_FC"][i].school == data["Top10_School_China_FC"][0][k])
          dataArr.push([Date.parse(data["Date_School_AU_China_FC"][i].date), parseInt(data["Date_School_AU_China_FC"][i]["Active Users"])])
      }
      dsacfc.push({
        name : data["Top10_School_China_FC"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_AU_China_FC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
          text: 'China FC Top10 School'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dsacfc
    });
  },
  // Week_School_AU_China_FC
  Week_School_AU_China_FC: function(data) {
    var wsacfc = [];
    for(var k in data["Top10_School_China_FC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Week_School_AU_China_FC"].length; i ++) {
        if(data["Week_School_AU_China_FC"][i].school == data["Top10_School_China_FC"][0][k])
          dataArr.push([Date.parse(data["Week_School_AU_China_FC"][i].date), parseInt(data["Week_School_AU_China_FC"][i]["Active Users"])])
      }
      wsacfc.push({
        name : data["Top10_School_China_FC"][0][k],
        data : dataArr
      })
    }
    $('#Week_School_AU_China_FC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
          text: 'China FC Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsacfc
    });
  },

  // Month_School_AU_China_FC
  Month_School_AU_China_FC: function(data) {
    var msacfc = [];
    for(var k in data["Top10_School_China_FC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Month_School_AU_China_FC"].length; i ++) {
        if(data["Month_School_AU_China_FC"][i].school == data["Top10_School_China_FC"][0][k])
          dataArr.push([Date.parse(data["Month_School_AU_China_FC"][i].date), parseInt(data["Month_School_AU_China_FC"][i]["Active Users"])])
      }
      msacfc.push({
        name : data["Top10_School_China_FC"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_AU_China_FC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
          text: 'China FC Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msacfc
    });
  },
    //Date_School_AU_China_SC
  Date_School_AU_China_SC: function(data) {
    var dsacsc = [];
    for(var k in data["Top10_School_China_SC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Date_School_AU_China_SC"].length; i ++) {
        if(data["Date_School_AU_China_SC"][i].school == data["Top10_School_China_SC"][0][k])
          dataArr.push([Date.parse(data["Date_School_AU_China_SC"][i].date), parseInt(data["Date_School_AU_China_SC"][i]["Active Users"])])
      }
      dsacsc.push({
        name : data["Top10_School_China_SC"][0][k],
        data : dataArr
      })
    }
    $('#Date_School_AU_China_SC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
          text: 'China SC Top10 School'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: dsacsc
    });
  },
  // Week_School_AU_China_SC
  Week_School_AU_China_SC: function(data) {
    var wsacsc = [];
    for(var k in data["Top10_School_China_SC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Week_School_AU_China_SC"].length; i ++) {
        if(data["Week_School_AU_China_SC"][i].school == data["Top10_School_China_SC"][0][k])
          dataArr.push([Date.parse(data["Week_School_AU_China_SC"][i].date), parseInt(data["Week_School_AU_China_SC"][i]["Active Users"])])
      }
      wsacsc.push({
        name : data["Top10_School_China_SC"][0][k],
        data : dataArr
      })
    }
    $('#Week_School_AU_China_SC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
          text: 'China SC Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsacsc
    });
  },

  // Month_School_AU_China_SC
  Month_School_AU_China_SC: function(data) {
    var msacsc = [];
    for(var k in data["Top10_School_China_SC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Month_School_AU_China_SC"].length; i ++) {
        if(data["Month_School_AU_China_SC"][i].school == data["Top10_School_China_SC"][0][k])
          dataArr.push([Date.parse(data["Month_School_AU_China_SC"][i].date), parseInt(data["Month_School_AU_China_SC"][i]["Active Users"])])
      }
      msacsc.push({
        name : data["Top10_School_China_SC"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_AU_China_SC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
          text: 'China SC Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 01, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msacsc
    });
  },

  // Week_School_AU_China_SC
  Week_School_AU_China_SC: function(data) {
    var wsacsc = [];
    for(var k in data["Top10_School_China_SC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Week_School_AU_China_SC"].length; i ++) {
        if(data["Week_School_AU_China_SC"][i].school == data["Top10_School_China_SC"][0][k])
          dataArr.push([Date.parse(data["Week_School_AU_China_SC"][i].date), parseInt(data["Week_School_AU_China_SC"][i]["Active Users"])])
      }
      wsacsc.push({
        name : data["Top10_School_China_SC"][0][k],
        data : dataArr
      })
    }
    $('#Week_School_AU_China_SC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
          text: 'China SC Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wsacsc
    });
  },

  // Month_School_AU_China_SC
  Month_School_AU_China_SC: function(data) {
    var msacsc = [];
    for(var k in data["Top10_School_China_SC"][0]) {
      var dataArr = [];
      for(var i = 0 ; i < data["Month_School_AU_China_SC"].length; i ++) {
        if(data["Month_School_AU_China_SC"][i].school == data["Top10_School_China_SC"][0][k])
          dataArr.push([Date.parse(data["Month_School_AU_China_SC"][i].date), parseInt(data["Month_School_AU_China_SC"][i]["Active Users"])])
      }
      msacsc.push({
        name : data["Top10_School_China_SC"][0][k],
        data : dataArr
      })
    }
    $('#Month_School_AU_China_SC').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
          text: 'China SC Top 10 Schools'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 01, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },
        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: msacsc
    });
  },

  //Date_Region_TU_Global
  Date_Region_TU_Global: function(data) {
    var drtu= [];
    for(var k in data["Top10_Region_Global"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_Region_TU_Global"].length; i ++) {
        if(data["Date_Region_TU_Global"][i].country == data["Top10_Region_Global"][0][k])
          dataArr.push([Date.parse(data["Date_Region_TU_Global"][i].date), parseInt(data["Date_Region_TU_Global"][i]["Total Users"])])
      }
      drtu.push({
        name : data["Top10_Region_Global"][0][k],
        data : dataArr
      })
    }
    $('#Date_Region_TU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily Total Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: drtu
    });
  },

  //Date_Region_NU_Global
  Date_Region_NU_Global: function(data) {
    var drnu= [];
    for(var k in data["Top10_Region_Global"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_Region_NU_Global"].length; i ++) {
        if(data["Date_Region_NU_Global"][i].country == data["Top10_Region_Global"][0][k])
          dataArr.push([Date.parse(data["Date_Region_NU_Global"][i].date), parseInt(data["Date_Region_NU_Global"][i]["New Users"])])
      }
      drnu.push({
        name : data["Top10_Region_Global"][0][k],
        data : dataArr
      })
    }
    $('#Date_Region_NU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily New Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: drnu
    });
  },

  //Week_Region_NU_Global
  Week_Region_NU_Global: function(data) {
    var wrnu= [];
    for(var k in data["Top10_Region_Global"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Week_Region_NU_Global"].length; i ++) {
        if(data["Week_Region_NU_Global"][i].country == data["Top10_Region_Global"][0][k])
          dataArr.push([Date.parse(data["Week_Region_NU_Global"][i].date), parseInt(data["Week_Region_NU_Global"][i]["New Users"])])
      }
      wrnu.push({
        name : data["Top10_Region_Global"][0][k],
        data : dataArr
      })
    }
    $('#Week_Region_NU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly New Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: wrnu
    });
  },

  //Month_Region_NU_Global
  Month_Region_NU_Global: function(data) {
    var mrnu= [];
    for(var k in data["Top10_Region_Global"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Month_Region_NU_Global"].length; i ++) {
        if(data["Month_Region_NU_Global"][i].country == data["Top10_Region_Global"][0][k])
          dataArr.push([Date.parse(data["Month_Region_NU_Global"][i].date), parseInt(data["Month_Region_NU_Global"][i]["New Users"])])
      }
      mrnu.push({
        name : data["Top10_Region_Global"][0][k],
        data : dataArr
      })
    }
    $('#Month_Region_NU_Global').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly New Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: mrnu
    });
  },

  //Date_Region_AU_Global_FX
  Date_Region_AU_Global_FX: function(data) {
    var drau= [];
    for(var k in data["Top10_Region_Global_FX"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Date_Region_AU_Global_FX"].length; i ++) {
        if(data["Date_Region_AU_Global_FX"][i].country == data["Top10_Region_Global_FX"][0][k])
          dataArr.push([Date.parse(data["Date_Region_AU_Global_FX"][i].date), parseInt(data["Date_Region_AU_Global_FX"][i]["Active Users"])])
      }
      drau.push({
        name : data["Top10_Region_Global_FX"][0][k],
        data : dataArr
      })
    }
    $('#Date_Region_AU_Global_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Daily Active Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: drau
    });
  },

  //Week_Region_AU_Global_FX
  Week_Region_AU_Global_FX: function(data) {
    var drau= [];
    for(var k in data["Top10_Region_Global_FX"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Week_Region_AU_Global_FX"].length; i ++) {
        if(data["Week_Region_AU_Global_FX"][i].country == data["Top10_Region_Global_FX"][0][k])
          dataArr.push([Date.parse(data["Week_Region_AU_Global_FX"][i].date), parseInt(data["Week_Region_AU_Global_FX"][i]["Active Users"])])
      }
      drau.push({
        name : data["Top10_Region_Global_FX"][0][k],
        data : dataArr
      })
    }
    $('#Week_Region_AU_Global_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Weekly Active Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: drau
    });
  },

  //Month_Region_AU_Global_FX
  Month_Region_AU_Global_FX: function(data) {
    var drau= [];
    for(var k in data["Top10_Region_Global_FX"][0]) {
      var dataArr = []
      for(var i = 0 ; i < data["Month_Region_AU_Global_FX"].length; i ++) {
        if(data["Month_Region_AU_Global_FX"][i].country == data["Top10_Region_Global_FX"][0][k])
          dataArr.push([Date.parse(data["Month_Region_AU_Global_FX"][i].date), parseInt(data["Month_Region_AU_Global_FX"][i]["Active Users"])])
      }
      drau.push({
        name : data["Top10_Region_Global_FX"][0][k],
        data : dataArr
      })
    }
    $('#Month_Region_AU_Global_FX').highcharts({
        chart: {
                backgroundColor: null,
                zoomType: 'x'
        },
        title: {
            text: 'Monthly Active Users'
        },
        subtitle: {
          text: 'Global Top10 Region'
        },
        xAxis: {
            type: 'datetime',
            dateTimeLabelFormats: { // don't display the dummy year
                month: '%b. %e',
                year: '%b'
            },
            title: {
                text: 'Date'
            },
            floor: Date.UTC(2015, 00, 01)
        },
        yAxis: {
            title: {
                text: 'Users'
            },
            min: 0
        },
        tooltip: {
            headerFormat: '<b>{series.name}</b><br>',
            pointFormat: '{point.x:%e. %b}: {point.y:.f} users'
        },

        plotOptions: {
            spline: {
                marker: {
                    enabled: false
                }
            }
        },
        credits: {
          enabled: false
        },
        series: drau
    });
  }




}
