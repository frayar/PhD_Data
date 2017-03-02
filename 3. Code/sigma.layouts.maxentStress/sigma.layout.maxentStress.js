;(function(undefined) {
  'use strict';

  if (typeof sigma === 'undefined')
    throw new Error('sigma is not declared');

  // Initialize package:
  sigma.utils.pkg('sigma.layouts.maxentStress');

  /**
   * Sigma Maxent Stress Majorization
   * ===============================
   *
   * Author: Frédéric RAYAR @ Université François-Rabelais of Tours
   * Version: 0.1
   * Date: 16-02-2016
   */


  /**
   * Settings Object
   * ------------------
   */
  var settings = {
    epsilon: 0.0001,		// 10e-4
    iterations: 250,		// maximum number of iterations
	alphaStep: 50, 			// steps per setting of alpha (after each step, apha *= 0.3 )
    q: 0
  };

  /**
   * Instance Object
   * ------------------
   */
  var _instance = {};


  /**
   * Event emitter Object
   * ------------------
   */
  var _eventEmitter = {};

  /**
   * Maxent Stress Majorization Object
   * ---------------------------
   */
  function MaxentStressMajorization() {
    var self = this;

	/**
     * Object initialization function.
     */
    this.init = function (sigInst, options) {
      options = options || {};
	
      // Properties
      this.sigInst = sigInst;
      this.config = sigma.utils.extend(options, settings);
      this.easing = options.easing;
      this.duration = options.duration;

      if (!sigma.plugins || typeof sigma.plugins.animate === 'undefined') {
        throw new Error('sigma.plugins.animate is not declared');
      }

      // State
      this.running = false;
    };


	/**
     * Compute Distance matrix between pairs of node.
     */
    this.computeDistancesMatrix = function () {

	  // Get nodes
	  var nodes = this.sigInst.graph.nodes(),
	      edges = this.sigInst.graph.edges(),
	      nodesCount = nodes.length,
		  i, j,
		  edge, src, tgt, pos;

	  // Initialize desired distances matrix
	  var distances = new Array();
	  for (i = 0; i < nodesCount; i++) {
		distances[i] = new Array();
		  for (j = 0; j < nodesCount; j++) 
			distances[i][j] = 0;
	  }		
	  
	  // Get weight of the edge between node1 and node2
	  for(i = 0; i<edges.length; i++) {
		edge = edges[i];
		src = edge.source; pos = src.lastIndexOf('.'); pos = Math.max(0,pos);  src = src.substring(pos + 1);
		tgt = edge.target; pos = tgt.lastIndexOf('.'); pos = Math.max(0,pos);  tgt = tgt.substring(pos + 1);
		distances[src][tgt] = distances[tgt][src] = edge.weight;
	  }

	  return distances;
	};	
	
    /**
     * Compute weights (w_ij)
     */
    this.computeWeight = function (distances) {
		 	  
	  // Get matrix dimension
	  var nodeCount = distances.length;
	  
	  // Initialize desired distances matrix
	  var weights = new Array();
	  for (var i = 0; i < nodeCount; i++) {
		weights[i] = new Array();
		  for (var j = 0; j < nodeCount; j++) 
			weights[i][j] = 0;
	  }
	  
	  // Compute desired distances
	  for (i = 0; i < nodeCount; i++)
		  for (j = i+1; j < nodeCount; j++) 
		    if (distances[i][j] != 0)
			  weights[i][j] = weights[j][i] = 1 / ( distances[i][j] * distances[i][j] );

	  return weights;
	};

    /**
     * Compute stress
     */
    this.computeStress = function (distances, weights) {
	
	  // Get nodes
	  var nodes = this.sigInst.graph.nodes(),
	      nodeCount = nodes.length,
	      stress = 0,
		  dist = 0,
		  n1, n2;
	  
	  // Compute stress
	  for (var i = 0; i < nodeCount; i++) {
	    n1 = nodes[i];
		for (var j = i+1; j < nodeCount; j++) 
		{
		  n2 = nodes[j];
		  // Compute euclidian distance
		  dist = (n1.sm_x - n2.sm_x) * (n1.sm_x - n2.sm_x) + (n1.sm_y - n2.sm_y) * (n1.sm_y - n2.sm_y);
		  dist = Math.sqrt( dist );
			
		  // Update stress
		  stress += weights[i][j] * (dist - distances[i][j]) * (dist - distances[i][j]);
		}
	  }

	  return stress;
	};	

    /**
     * Single layout iteration.
     */
    this.atomicGo = function (distances, weights) {
      if (!this.running || this.iterCount < 1) return false;

      var nodes = this.sigInst.graph.nodes(),
	      nodesCount = nodes.length,
          n1,
          n2,
          dist = 0.0,
		  newX = 0.0,
		  newY = 0.0,
		  totalWeight = 0.0,
		  alpha = 0.0,
		  force = 0.0;
		  
	  this.iterCount--;
      this.running = (this.iterCount > 0);
	  
	  // Compute alpha value
	  alpha = Math.pow( 0.3, Math.floor((this.config.iterations-this.iterCount)/this.config.alphaStep) );
	  alpha = this.config.q >= 0 ? alpha : -alpha;

      for (var i = 0; i < nodesCount; i++) {
        n1 = nodes[i];
		newX = 0.0;
		newY = 0.0;
		totalWeight = 0.0;

        for (var j = 0; j < nodesCount; j++) {
          n2 = nodes[j];
		  
		  // Do not consider
		  if (n1.id == n2.id) continue;
		  
		  // Compute euclidian distance
		  dist = (n1.sm_x - n2.sm_x) * (n1.sm_x - n2.sm_x) + (n1.sm_y - n2.sm_y) * (n1.sm_y - n2.sm_y);
		  dist = Math.sqrt( dist );
		  
		  // Check if an edge ij exist in the graph
		  if ( distances[i][j] != 0 ) {

			// Compute n2 contribution to n1 on the X-axis
			newX += weights[i][j] * ( ( n2.sm_x + ( distances[i][j] * (n1.sm_x - n2.sm_x) / dist ) ) );
		  
			// Compute n2 contribution to n1 on the Y-axis
			newY += weights[i][j] * ( ( n2.sm_y + ( distances[i][j] * (n1.sm_y - n2.sm_y) / dist ) ) );
		  
			// Update total weight
			totalWeight += weights[i][j]
			
		  } else { 	// no edge, consider attraction or repulsion force depending whether q < 0 or q >=0
		  
			// compute force 
			force = Math.pow(dist,this.config.q + 2)
			
		    // Compute n2 contribution to n1 on the X-axis
			newX += alpha * (n1.sm_x - n2.sm_x) / force ;
			
		    // Compute n2 contribution to n1 on the Y-axis
			newY += alpha * (n1.sm_y - n2.sm_y) / force ;			
		  }
		}

		// Update n1 coordinates
		n1.sm_x = newX / totalWeight;
		n1.sm_y = newY / totalWeight;
		
		// Debug
		/*
		console.log("newX = " + newX + " - newY = " + newY + " - totalWeight = " + totalWeight );
		console.log("Post -Iteration #" + this.iterCount + " = (" + n1.sm_x  + "," + n1.sm_y + ")");
		*/

	  }

      return this.running;
    };

	/**
     * Check if the stop criterion is met
     */
    this.finished = function () {
	  return ( this.previousStress - this.currentStress < this.previousStress * this.config.epsilon )
	};
	
    /**
     * Layout iterations.
     */
    this.go = function (distances, weights) {
	
	  this.iterCount = this.config.iterations;

	  // Calculate stress for the first time
	  this.currentStress = this.computeStress(distances, weights);
	    
	  // Do-while loop
	  do {

		// Process an iteration
        this.atomicGo(distances, weights);
		
		// Update stress value
		this.previousStress = this.currentStress;
		this.currentStress = this.computeStress(distances, weights);
		
	  //} while ( !this.finished() && this.running);
	  } while (this.running);
	  
	  // Debug
	  /*
	  var numberIter = this.config.iterations - this.iterCount;
	  var achieved_epsilon = (this.previousStress - this.currentStress) / this.previousStress;
	  console.log( 'Processed in ' +  numberIter + ' iterations. Achieved epsilon = ' + achieved_epsilon );
	  console.log('this.previousStress = ' + this.previousStress);
	  console.log('this.currentStress = ' + this.currentStress);
	  */
	  
	  // Debug
	  /*
	  var nodes = this.sigInst.graph.nodes();
	  for (var i = 0; i < nodes.length; i++)
          console.log("Node #" + i + " = (" + nodes[i].sm_x  + "," + nodes[i].sm_y + ")");
	  */

				
	  // Stop layout computation
      this.stop();
    };
	
    /**
     * Layout computation.
     */
    this.start = function() {
      if (this.running) return;

      var nodes = this.sigInst.graph.nodes();

      this.running = true;

      // Init nodes stress majorization coordintates
      for (var i = 0; i < nodes.length; i++) {
        nodes[i].sm_x = nodes[i].x;
        nodes[i].sm_y = nodes[i].y;
      }
	  
	  // Compute shortest path as desired distances
	  var distances = this.computeDistancesMatrix();
	  
	  // Debug
	  /*
	  for (var i = 0; i < distances.length; i++) {
	    var msg = "";
		for (var j = 0; j < distances.length; j++) 
		  msg += distances[i][j] + " ";
		console.log(msg);
	  } 
	  */

	  // Compute weights
	  var weights = this.computeWeight(distances);
	  
      _eventEmitter[self.sigInst.id].dispatchEvent('start');

	  // Go
      this.go(distances, weights);
    };

	/**
     * Finish layout computation.
     */
    this.stop = function() {
      var nodes = this.sigInst.graph.nodes();

      this.running = false;

      if (this.easing) {
        _eventEmitter[self.sigInst.id].dispatchEvent('interpolate');
        sigma.plugins.animate(
          self.sigInst,
          {
            x: 'sm_x',
            y: 'sm_y'
          },
          {
            easing: self.easing,
            onComplete: function() {
              self.sigInst.refresh();
              for (var i = 0; i < nodes.length; i++) {
                nodes[i].sm_x = null;
                nodes[i].sm_y = null;
              }
              _eventEmitter[self.sigInst.id].dispatchEvent('stop');
            },
            duration: self.duration
          }
        );
      }
      else {
        // Apply changes
        for (var i = 0; i < nodes.length; i++) {
          nodes[i].x = nodes[i].sm_x;
          nodes[i].y = nodes[i].sm_y;
        }

        this.sigInst.refresh();

        for (var i = 0; i < nodes.length; i++) {
          nodes[i].sm_x = null;
          nodes[i].sm_y = null;
        }
        _eventEmitter[self.sigInst.id].dispatchEvent('stop');
      }
    };

	/**
     * Kill layout computation.
     */
    this.kill = function() {
      this.sigInst = null;
      this.config = null;
      this.easing = null;
    };
  };




  /**
   * Interface
   * ----------
   */

  /**
   * Configure the layout algorithm.

   * Recognized options:
   * **********************
   * Here is the exhaustive list of every accepted parameters in the settings
   * object:
   *
   *   {?boolean}           autoArea   If `true`, area will be computed as N².
   *   {?number}            area       The area of the graph.
   *   {?number}            gravity    This force attracts all nodes to the
   *                                   center to avoid dispersion of
   *                                   disconnected components.
   *   {?number}            speed      A greater value increases the
   *                                   convergence speed at the cost of precision loss.
   *   {?number}            iterations The number of iterations to perform
   *                                   before the layout completes.
   *   {?(function|string)} easing     Either the name of an easing in the
   *                                   sigma.utils.easings package or a
   *                                   function. If not specified, the
   *                                   quadraticInOut easing from this package
   *                                   will be used instead.
   *   {?number}            duration   The duration of the animation. If not
   *                                   specified, the "animationsTime" setting
   *                                   value of the sigma instance will be used
   *                                   instead.
   *
   *
   * @param  {sigma}   sigInst The related sigma instance.
   * @param  {object} config  The optional configuration object.
   *
   * @return {sigma.classes.dispatcher} Returns an event emitter.
   */
  sigma.layouts.maxentStress.configure = function(sigInst, config) {
    if (!sigInst) throw new Error('Missing argument: "sigInst"');
    if (!config) throw new Error('Missing argument: "config"');

    // Create instance if undefined
    if (!_instance[sigInst.id]) {
      _instance[sigInst.id] = new MaxentStressMajorization();

      _eventEmitter[sigInst.id] = {};
      sigma.classes.dispatcher.extend(_eventEmitter[sigInst.id]);

      // Binding on kill to clear the references
      sigInst.bind('kill', function() {
        _instance[sigInst.id].kill();
        _instance[sigInst.id] = null;
        _eventEmitter[sigInst.id] = null;
      });
    }

    _instance[sigInst.id].init(sigInst, config);

    return _eventEmitter[sigInst.id];
  };

  /**
   * Start the layout algorithm. It will use the existing configuration if no
   * new configuration is passed.

   * Recognized options:
   * **********************
   * Here is the exhaustive list of every accepted parameters in the settings
   * object:
   *
   *   {?boolean}           autoArea   If `true`, area will be computed as N².
   *   {?number}            area       The area of the graph.
   *   {?number}            gravity    This force attracts all nodes to the
   *                                   center to avoid dispersion of
   *                                   disconnected components.
   *   {?number}            speed      A greater value increases the
   *                                   convergence speed at the cost of precision loss.
   *   {?number}            iterations The number of iterations to perform
   *                                   before the layout completes.
   *   {?(function|string)} easing     Either the name of an easing in the
   *                                   sigma.utils.easings package or a
   *                                   function. If not specified, the
   *                                   quadraticInOut easing from this package
   *                                   will be used instead.
   *   {?number}            duration   The duration of the animation. If not
   *                                   specified, the "animationsTime" setting
   *                                   value of the sigma instance will be used
   *                                   instead.
   *
   *
   * @param  {sigma}   sigInst The related sigma instance.
   * @param  {?object} config  The optional configuration object.
   *
   * @return {sigma.classes.dispatcher} Returns an event emitter.
   */
  sigma.layouts.maxentStress.start = function(sigInst, config) {
    if (!sigInst) throw new Error('Missing argument: "sigInst"');

    if (config) {
      this.configure(sigInst, config);
    }

    _instance[sigInst.id].start();

    return _eventEmitter[sigInst.id];
  };

  /**
   * Returns true if the layout has started and is not completed.
   *
   * @param  {sigma}   sigInst The related sigma instance.
   *
   * @return {boolean}
   */
  sigma.layouts.maxentStress.isRunning = function(sigInst) {
    if (!sigInst) throw new Error('Missing argument: "sigInst"');

    return !!_instance[sigInst.id] && _instance[sigInst.id].running;
  };

  /**
   * Returns the number of iterations done divided by the total number of
   * iterations to perform.
   *
   * @param  {sigma}   sigInst The related sigma instance.
   *
   * @return {number} A value between 0 and 1.
   */
  sigma.layouts.maxentStress.progress = function(sigInst) {
    if (!sigInst) throw new Error('Missing argument: "sigInst"');

    return (_instance[sigInst.id].config.iterations - _instance[sigInst.id].iterCount) /
      _instance[sigInst.id].config.iterations;
  };
  
}).call(this);
