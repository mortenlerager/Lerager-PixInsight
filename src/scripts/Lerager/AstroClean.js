#engine v8
#feature-id Lerager > AstroClean
#feature-info AstroClean v2.8.7. Correct a linear image with a synthetic flat estimated from its starless counterpart.

#include <pjsr/StdButton.jsh>
#include <pjsr/StdIcon.jsh>

"use strict";

const AC_DEBUG_MASK = false; // Development build: false disables diagnostic images.
function acSpinner(label) {
   let phase=0;
   const timer=new Timer;timer.interval=0.12;timer.periodic=true;
   timer.onTimeout=()=>{label.text=["|","/","—","\\"][phase++%4]+"  "+label.text.replace(/^.  /, "");};
   timer.start();return timer;
}
const AC_FLAT_FACTOR = 4; // Internal only: 4=25%, 2=50%, 1=100%.

function acSelect(values, k, length = values.length) {
   let lo = 0, hi = length - 1;
   while (lo < hi) {
      const pivot = values[(lo + hi) >> 1];
      let i = lo, j = hi;
      while (i <= j) {
         while (values[i] < pivot) ++i;
         while (values[j] > pivot) --j;
         if (i <= j) {
            const tmp = values[i];
            values[i++] = values[j];
            values[j--] = tmp;
         }
      }
      if (k <= j) hi = j;
      else if (k >= i) lo = i;
      else break;
   }
   return values[k];
}

function acMedian(values) {
   const copy = values.slice();
   const k = copy.length >> 1;
   const upper = acSelect(copy, k);
   return copy.length % 2 ? upper : (upper + acSelect(copy, k-1))/2;
}

function acStats(values) {
   const median = acMedian(values);
   const dev = values.map(value => Math.abs(value - median));
   return [median, Math.max(1e-12, 1.4826 * acMedian(dev))];
}

function acReflect(i, n) {
   if (n === 1) return 0;
   while (i < 0 || i >= n)
      i = i < 0 ? -1-i : 2*n-1-i;
   return i;
}

function acFilter(src, w, h, r, median, tick) {
   const size = 2*r + 1;
   const out = new Float32Array(src.length);
   const kernel = new Float32Array(size*size);
   const xs = new Int32Array(w*size);
   for (let x = 0; x < w; ++x)
      for (let k = 0; k < size; ++k)
         xs[x*size+k] = acReflect(x+k-r,w);
   for (let y = 0; y < h; ++y) {
      if (y % 32 === 0) tick();
      const rows = [];
      for (let k = -r; k <= r; ++k) rows.push(acReflect(y+k,h)*w);
      for (let x = 0; x < w; ++x) {
         let i = 0, sum = 0;
         for (const row of rows)
            for (let k = 0; k < size; ++k) {
               const value = src[row + xs[x*size+k]];
               if (median) kernel[i++] = value;
               else sum += value;
            }
         out[y*w+x] = median ? acSelect(kernel,kernel.length >> 1) : sum/kernel.length;
      }
   }
   return out;
}

function acHalf(src, w, h) {
   const width = w >> 1, height = h >> 1;
   const out = new Float32Array(width*height);
   for (let y = 0; y < height; ++y)
      for (let x = 0; x < width; ++x) {
         const i = 2*(y*w+x);
         out[y*width+x] = (src[i]+src[i+1]+src[i+w]+src[i+w+1])/4;
      }
   return out;
}

function acUpsample(src, w, h, width, height, binary = false) {
   if(w===width && h===height) return src;
   const out = binary ? new Uint8Array(width*height) : new Float32Array(width*height);
   const sx = (w-1)/(width-1), sy = (h-1)/(height-1);
   for (let y = 0; y < height; ++y) {
      const yy = y*sy, iy = Math.floor(yy), fy = yy-iy;
      const row = iy*w, next = Math.min(iy+1,h-1)*w;
      for (let x = 0; x < width; ++x) {
         const xx = x*sx, ix = Math.floor(xx), fx = xx-ix;
         if (binary) out[y*width+x] = src[Math.round(yy)*w+Math.round(xx)];
         else {
            const j = Math.min(ix+1,w-1);
            const a = src[row+ix]*(1-fx)+src[row+j]*fx;
            const b = src[next+ix]*(1-fx)+src[next+j]*fx;
            out[y*width+x] = a*(1-fy)+b*fy;
         }
      }
   }
   return out;
}

function acPyramid(src, w, h, levels, tick) {
   const result = [src];
   let width = w, height = h;
   for (let i = 0; i < levels; ++i) {
      src = acHalf(acFilter(src,width,height,2,true,tick),width,height);
      width >>= 1;
      height >>= 1;
      result.push(acUpsample(src,width,height,w,h));
   }
   return result;
}

function acNeighbors(p,w,h,diagonal,visit) {
   const x = p%w, y = Math.floor(p/w);
   if (x > 0) visit(p-1);
   if (x+1 < w) visit(p+1);
   if (y > 0) visit(p-w);
   if (y+1 < h) visit(p+w);
   if (diagonal) {
      if (y > 0 && x > 0) visit(p-w-1);
      if (y > 0 && x+1 < w) visit(p-w+1);
      if (y+1 < h && x > 0) visit(p+w-1);
      if (y+1 < h && x+1 < w) visit(p+w+1);
   }
}

function acComponents(mask,w,h,consume) {
   const todo = mask.slice();
   const queue = new Int32Array(mask.length);
   for (let p = 0; p < todo.length; ++p) {
      if (!todo[p]) continue;
      let size = 1;
      queue[0] = p;
      todo[p] = 0;
      for (let i = 0; i < size; ++i)
         acNeighbors(queue[i],w,h,true,q => {
            if (todo[q]) {
               todo[q] = 0;
               queue[size++] = q;
            }
         });
      consume(queue.subarray(0,size));
   }
}

function acFlood(seeds,support,w,h) {
   const flooded = seeds.slice(), queue = new Int32Array(seeds.length);
   let size = 0;
   for (let i = 0; i < seeds.length; ++i)
      if (seeds[i]) queue[size++] = i;
   for (let i = 0; i < size; ++i)
      acNeighbors(queue[i],w,h,false,q => {
         if (support[q] && !flooded[q]) {
            flooded[q] = 1;
            queue[size++] = q;
         }
      });
   return flooded;
}

function acDilate(mask,w,h,r) {
   const result = mask.slice();
   for (let p = 0; p < mask.length; ++p) {
      if (!mask[p]) continue;
      const x = p%w, y = Math.floor(p/w);
      if (x > 0 && x+1 < w && y > 0 && y+1 < h &&
          mask[p-1] && mask[p+1] && mask[p-w] && mask[p+w]) continue;
      for (let dy = Math.max(-r,-y); dy <= Math.min(r,h-1-y); ++dy) {
         const dx = Math.floor(Math.sqrt(r*r-dy*dy));
         const center = p+dy*w;
         result.fill(1,center+Math.max(-dx,-x),center+Math.min(dx,w-1-x)+1);
      }
   }
   return result;
}

// Reject complete dark candidates and dark peripheral lobes, not individual pixels.
function acRejectDarkObjects(coarse,core,small,w,h,tick) {
   const local=acFilter(coarse,w,h,8,true,tick);
   const delta=coarse.map((v,i)=>v-local[i]);
   const [dm,ds]=acStats(delta);
   const dark=new Uint8Array(coarse.length);
   acComponents(Uint8Array.from(delta,v=>v<Math.min(0,dm-3*ds)?1:0),w,h,component=>{
      if(component.length<3 || component.length>500) return;
      // Preserve dark structures embedded in the detected bright object core.
      if(component.some(i=>core[i])) return;
      for(const i of component) dark[i]=1;
   });
   
   const accepted=new Uint8Array(small.length);
   acComponents(small,w,h,component=>{
      const footprint=new Set(), ring=new Set();
      for(const p of component){
         const x=p%w,y=Math.floor(p/w);
         for(let dy=-12;dy<=12;++dy) for(let dx=-12;dx<=12;++dx){
            const d=dx*dx+dy*dy;
            if(d>144 || x+dx<0 || x+dx>=w || y+dy<0 || y+dy>=h) continue;
            const i=p+dy*w+dx;
            if(d<=16) footprint.add(i); else if(d>=64) ring.add(i);
         }
      }
      const sky=Array.from(ring).filter(i=>!footprint.has(i)&&!core[i]).map(i=>coarse[i]);
      if(!sky.length) return;
      const baseline=acMedian(sky);
      const values=Array.from(footprint).map(i=>coarse[i]);
      let negative=0,positive=0;
      for(const i of footprint){
         if(delta[i]<Math.min(0,dm-2.5*ds)) ++negative;
         if(delta[i]>Math.max(0,dm+2.5*ds)) ++positive;
      }
      if(acMedian(values)>=baseline && positive>negative && !Array.from(footprint).some(i=>dark[i]))
         for(const i of component) accepted[i]=1;
   });
   const expanded=acDilate(core,w,h,8);
   const compact=acDilate(accepted,w,h,4);
   // Preserve continuous masks: do not carve dark holes after dilation.
   return Uint8Array.from(expanded,(v,i)=>v || compact[i] ? 1:0);
}

// Diagnostic background: lower-quartile tile estimates, smoothed on a sparse grid.
function acDiagnosticBackground(src,w,h,tick) {
   const nx=Math.max(2,Math.min(8,Math.floor(w/32)));
   const ny=Math.max(2,Math.min(8,Math.floor(h/32)));
   const grid=new Float32Array(nx*ny);
   for(let gy=0;gy<ny;++gy) {
      tick();
      for(let gx=0;gx<nx;++gx) {
         const values=[];
         const x0=Math.floor(gx*w/nx),x1=Math.floor((gx+1)*w/nx);
         const y0=Math.floor(gy*h/ny),y1=Math.floor((gy+1)*h/ny);
         for(let y=y0;y<y1;++y) for(let x=x0;x<x1;++x) values.push(src[y*w+x]);
         grid[gy*nx+gx]=acSelect(values,Math.floor(values.length*0.25));
      }
   }
   // Suppress bright contaminated tiles using a neighborhood lower quartile.
   const clean=grid.map((v,i)=>{
      const x=i%nx,y=Math.floor(i/nx),a=[];
      for(let yy=Math.max(0,y-1);yy<=Math.min(ny-1,y+1);++yy)
         for(let xx=Math.max(0,x-1);xx<=Math.min(nx-1,x+1);++xx) a.push(grid[yy*nx+xx]);
      return Math.min(v,acSelect(a,Math.floor(a.length*0.25)));
   });
   return acUpsample(acFilter(clean,nx,ny,1,false,tick),nx,ny,w,h);
}

function acMask(fine,fw,fh,tick,sensitivity=50,debug=AC_DEBUG_MASK) {
   const high=Math.max(0,(sensitivity-50)/50);
   const k=Math.pow(2,(50-sensitivity)/50-2*high);
   const w = fw >> 1, h = fh >> 1;
   const coarse = acHalf(fine,fw,fh);
   const levels = acPyramid(coarse,w,h,6,tick);
   const extended = levels[2].map((v,i) => v-levels[5][i]);
   const broad = levels[3].map((v,i) => v-levels[6][i]);
   const [em,es] = acStats(extended), [bm,bs] = acStats(broad);
   // Require corroboration at the smaller scale before broad-scale growth.
   // Broad residuals alone can connect dust-shadow halos to a bright galaxy.
   const support = Uint8Array.from(extended,(v,i) =>
      v > em+3.5*es*k || (v > em+2.5*es*k && broad[i] > bm+2.8*bs*k) ||
      (high>0 && v > em+(2.5-2.5*high)*es*k && broad[i] > bm+(5-3*high)*bs*k));
   const seeds = new Uint8Array(coarse.length);
   acComponents(Uint8Array.from(extended,v => v > em+10*es*k),w,h,component => {
      if (component.length >= 50 && component.some(i => extended[i] > 18*es*k))
         for (const i of component) seeds[i] = 1;
   });
   const alt=acDiagnosticBackground(coarse,w,h,tick);
   const signal=levels[2].map((v,i)=>v-alt[i]);
   const [am,as]=acStats(signal);
   const altSupport=Uint8Array.from(signal,v=>v>am+3.5*as*k);
   const grown=acFlood(seeds,altSupport,w,h);
   const open = grown.map(v => 1-v);
   const boundary = open.map((v,i) => v && (i < w || i >= w*(h-1) || i%w === 0 || i%w === w-1));
   const closed = acFlood(boundary,open,w,h).map(v => 1-v);
   const large = acDilate(closed,w,h,8);

   const smooth = acFilter(fine,fw,fh,1,false,tick);
   const fineMedian = acFilter(smooth,fw,fh,5,true,tick);
   const narrow = acHalf(smooth.map((v,i) => v-fineMedian[i]),fw,fh);
   const compact = acFilter(coarse,w,h,1,false,tick).map((v,i) => v-levels[2][i]);
   const [cm,cs] = acStats(compact), [nm,ns] = acStats(narrow);
   const candidates = Uint8Array.from(compact,(v,i) => {
      const x = i%w, y = Math.floor(i/w);
      return !large[i] && x >= 20 && x < w-20 && y >= 20 && y < h-20 && v > cm+4*cs;
   });
   const small = new Uint8Array(coarse.length);
   acComponents(candidates,w,h,component => {
      if (component.length >= 4 && component.length <= 250 &&
          component.some(i => compact[i] > cm+5.5*cs) &&
          component.some(i => narrow[i] >= nm+5*ns))
         for (const i of component) small[i] = 1;
   });
   const finalMask = acRejectDarkObjects(coarse,closed,small,w,h,tick);
   // Diagnostic only: export each stage at the original mask grid resolution.
   // Signed residuals: zero relative to the robust median is mid-gray.
   // Use the same sigma-relative display range for both residual displays, without changing detection.
   const residualDisplay=(data,median,sigma)=>Float32Array.from(data,v=>Math.max(0,Math.min(1,0.5+(v-median)/(40*sigma))));
   const stages = [];
   if (debug) {
      const [bgm,bgs]=acStats(levels[5]);
      stages.push(["01_OriginalFlood",acFlood(seeds,support,w,h)],
                  ["02_AlternativeBackground",residualDisplay(alt,bgm,bgs)],
                  ["03_AlternativeSignal",residualDisplay(signal,am,as)],
                  ["04_AlternativeSupport",altSupport],
                  ["05_AlternativeFlood",grown]);

   }
   if (debug) {
      // Replace previous diagnostic windows only; never touch source or result views.
      for (const win of ImageWindow.windows.slice())
         if (/^AC_Mask_(?:S\d+_)?\d{2}_(?:Seeds|Flood|Closed|Dilated|SmallCandidates|Final|Background_L5|Extended_L2_L5|Broad_L3_L6|Background_L6|FloodSupport|AlternativeBackground|AlternativeSignal|AlternativeSupport|AlternativeFlood|OriginalFlood)(?:_\d+)?$/.test(win.mainView.id))
            win.forceClose();
   }
   if (debug) for (const [name, data] of stages) {
      tick();
      const window = acCreateOutput("AC_Mask_S"+sensitivity+"_"+name,w,h,1);
      try {
         acEditOutput(window,image=>{
            for(let y=0;y<h;y+=64) {
               const rows=Math.min(64,h-y);
               acWriteRows(image,y,rows,0,Float32Array.from(data.subarray(y*w,(y+rows)*w)));
            }
         });
         window.show();
      } catch(e) {window.forceClose();throw e;}
   }
   return finalMask;
}

function acReconstruct(plane,mask,w,h,tick) {
   const filtered = acFilter(plane,w,h,1,true,tick);
   const output = Float64Array.from(plane), indices = [];
   const neighbor = [];
   let sum = 0, count = 0;
   for (let p = 0; p < mask.length; ++p) {
      if (!mask[p]) continue;
      indices.push(p);
      const adjacent = [];
      acNeighbors(p,w,h,false,q => {
         adjacent.push(q);
         if (!mask[q]) {
            output[q] = filtered[q];
            sum += filtered[q];
            ++count;
         }
      });
      neighbor.push(adjacent);
   }
   if (indices.length === 0) return plane.slice();
   if (!count) throw new Error("Object mask covers the entire image.");
   const n = indices.length;
   const positions = new Int32Array(mask.length).fill(-1);
   for (let i = 0; i < n; ++i) positions[indices[i]] = i;
   const adjacency = new Int32Array(n*4).fill(-1);
   const degree = new Float64Array(n);
   const rhs = new Float64Array(n);
   for (let i = 0; i < n; ++i) {
      degree[i] = neighbor[i].length;
      for (let k = 0; k < neighbor[i].length; ++k) {
         const q = neighbor[i][k];
         adjacency[4*i+k] = positions[q];
         if (!mask[q]) rhs[i] += filtered[q];
      }
   }
   function multiply(x,result) {
      for (let i = 0; i < n; ++i) {
         let value = degree[i]*x[i];
         for (let k = 0; k < 4; ++k) {
            const j = adjacency[4*i+k];
            if (j >= 0) value -= x[j];
         }
         result[i] = value;
      }
   }
   const x = new Float64Array(n).fill(sum/count);
   const residual = new Float64Array(n), direction = new Float64Array(n);
   const product = new Float64Array(n);
   multiply(x,product);
   let rz = 0;
   for (let i = 0; i < n; ++i) {
      residual[i] = rhs[i]-product[i];
      direction[i] = residual[i]/degree[i];
      rz += residual[i]*direction[i];
   }
   const tolerance = Math.max(Math.abs(sum/count)*1e-8,1e-12);
   let converged = false;
   for (let iteration = 0; iteration < 4000; ++iteration) {
      if (iteration % 10 === 0) tick("Solving masked area — iteration " + iteration);
      let maximum = 0;
      for (let i = 0; i < n; ++i) maximum = Math.max(maximum,Math.abs(residual[i]));
      if (maximum < tolerance) { converged = true; break; }
      multiply(direction,product);
      let curvature = 0;
      for (let i = 0; i < n; ++i) curvature += direction[i]*product[i];
      if (!(curvature > 0)) throw new Error("Background solver encountered an invalid system.");
      const alpha = rz/curvature;
      let nextRz = 0;
      for (let i = 0; i < n; ++i) {
         x[i] += alpha*direction[i];
         residual[i] -= alpha*product[i];
         nextRz += residual[i]*residual[i]/degree[i];
      }
      const beta = nextRz/rz;
      for (let i = 0; i < n; ++i)
         direction[i] = residual[i]/degree[i]+beta*direction[i];
      rz = nextRz;
   }
   if (!converged) throw new Error("Background reconstruction did not converge.");
   for (let i = 0; i < n; ++i) output[indices[i]] = x[i];
   return plane.map((v,i) => mask[i] ? output[i] : v);
}

// Independent local background interpolation. No BlemishBlaster code is embedded.
// Sources always come from the original starless plane, never from filled pixels.
function acLocalFill(plane,mask,w,h,regions,tick,sourceExclusion=mask) {
   const labels = new Int32Array(mask.length).fill(-1), objects=[];
   acComponents(mask,w,h,indices => {
      const id=objects.length;
      let x0=w,y0=h,x1=0,y1=0;
      for(const p of indices) {
         labels[p]=id;
         const x=p%w,y=Math.floor(p/w);
         x0=Math.min(x0,x); x1=Math.max(x1,x);
         y0=Math.min(y0,y); y1=Math.max(y1,y);
      }
      objects.push({x0,y0,x1,y1,regions:[]});
   });
   if(!objects.length) return plane.slice();
   // Associate each manual source with the nearest object boundary.
   for(const region of (regions || [])) {
      let best=-1,dBest=Infinity;
      for(let id=0;id<objects.length;++id) {
         const o=objects[id]; let d=Infinity;
         for(let p=0;p<region.length;++p) if(region[p] && !mask[p]) {
            const x=p%w,y=Math.floor(p/w);
            const dx=Math.max(o.x0-x,0,x-o.x1),dy=Math.max(o.y0-y,0,y-o.y1);
            d=Math.min(d,dx*dx+dy*dy);
         }
         if(d<dBest) {dBest=d;best=id;}
      }
      if(best>=0) objects[best].regions.push(region);
   }
   // Keep automatic donor tiles away from the mask's potentially contaminated edge.
   const excluded=sourceExclusion.slice();
   let edge=sourceExclusion.slice();
   for(let k=0;k<8;++k) {
      const next=new Uint8Array(mask.length);
      for(let p=0;p<edge.length;++p) if(edge[p])
         acNeighbors(p,w,h,true,q=>{if(!excluded[q]){excluded[q]=1;next[q]=1;}});
      edge=next;
   }
   function donors(allowed,manual) {
      const points=[],step=12;
      for(let y0=0;y0<h;y0+=step) {
         if(y0%48===0) tick('Measuring clean background');
         for(let x0=0;x0<w;x0+=step) {
            const values=[];let sx=0,sy=0;
            for(let y=y0;y<Math.min(h,y0+step);++y)
               for(let x=x0;x<Math.min(w,x0+step);++x) {
                  const p=y*w+x;
                  if(!sourceExclusion[p] && (manual ? allowed[p] : !excluded[p]) && isFinite(plane[p])) {
                     values.push(plane[p]);sx+=x;sy+=y;
                  }
               }
            if(values.length>=9) points.push([sx/values.length,sy/values.length,acMedian(values)]);
         }
      }
      return points;
   }
   function tree(points,depth=0) {
      if(!points.length) return null;
      const axis=depth%2;points.sort((a,b)=>a[axis]-b[axis]);
      const m=points.length>>1;
      return {p:points[m],axis,left:tree(points.slice(0,m),depth+1),right:tree(points.slice(m+1),depth+1)};
   }
   function estimate(root,x,y) {
      const near=[];
      function visit(node) {
         if(!node) return;
         const p=node.p,dx=x-p[0],dy=y-p[1],d=dx*dx+dy*dy;
         near.push([d,p[2]]);near.sort((a,b)=>a[0]-b[0]);if(near.length>12) near.pop();
         const delta=node.axis===0?dx:dy;
         visit(delta<0?node.left:node.right);
         if(near.length<12 || delta*delta<near[near.length-1][0]) visit(delta<0?node.right:node.left);
      }
      visit(root);
      let total=0,weight=0;
      for(const s of near){const a=1/(s[0]+144);total+=a*s[1];weight+=a;}
      return total/weight;
   }
   let automatic=null;
   const trees=objects.map(o=>{
      if(o.regions.length) {
         const allowed=new Uint8Array(mask.length);
         for(const region of o.regions) for(let p=0;p<mask.length;++p) allowed[p]|=region[p];
         const samples=donors(allowed,true);
         if(!samples.length) throw new Error('Background source is too small. Enlarge its lasso.');
         return tree(samples);
      }
      if(!automatic) {
         const samples=donors(null,false);
         if(!samples.length) throw new Error('Not enough clear background. Mark a source outside the object mask.');
         automatic=tree(samples);
      }
      return automatic;
   });
   const result=plane.slice();
   // Small spatial blocks limit working memory. There is no copied texture and
   // no feather containing the original galaxy light inside the object mask.
   for(let by=0;by<h;by+=16) {
      tick('Filling object areas — '+Math.round(100*by/h)+'%');
      for(let bx=0;bx<w;bx+=16)
         for(let y=by;y<Math.min(by+16,h);++y)
            for(let x=bx;x<Math.min(bx+16,w);++x) {
               const p=y*w+x;
               if(mask[p]) result[p]=estimate(trees[labels[p]],x,y);
            }
   }
   return result;
}

function acAlternativeFill(plane,mask,w,h,regions,tick,sourceExclusion=mask) {
   const labels = new Int32Array(mask.length).fill(-1), objects=[];
   acComponents(mask,w,h,indices => {
      const id=objects.length;
      let x0=w,y0=h,x1=0,y1=0;
      for(const p of indices) {
         labels[p]=id;
         const x=p%w,y=Math.floor(p/w);
         x0=Math.min(x0,x); x1=Math.max(x1,x);
         y0=Math.min(y0,y); y1=Math.max(y1,y);
      }
      objects.push({x0,y0,x1,y1,regions:[]});
   });
   if(!objects.length) return plane.slice();
   // Associate each manual source with the nearest object boundary.
   for(const region of (regions || [])) {
      let best=-1,dBest=Infinity;
      for(let id=0;id<objects.length;++id) {
         const o=objects[id]; let d=Infinity;
         for(let p=0;p<region.length;++p) if(region[p] && !mask[p]) {
            const x=p%w,y=Math.floor(p/w);
            const dx=Math.max(o.x0-x,0,x-o.x1),dy=Math.max(o.y0-y,0,y-o.y1);
            d=Math.min(d,dx*dx+dy*dy);
         }
         if(d<dBest) {dBest=d;best=id;}
      }
      if(best>=0) objects[best].regions.push(region);
   }
   // Keep automatic donor tiles away from the mask's potentially contaminated edge.
   const excluded=sourceExclusion.slice();
   let edge=sourceExclusion.slice();
   for(let k=0;k<8;++k) {
      const next=new Uint8Array(mask.length);
      for(let p=0;p<edge.length;++p) if(edge[p])
         acNeighbors(p,w,h,true,q=>{if(!excluded[q]){excluded[q]=1;next[q]=1;}});
      edge=next;
   }
   function donors(allowed,manual) {
      const points=[],step=12;
      for(let y0=0;y0<h;y0+=step) {
         if(y0%48===0) tick('Measuring clean background');
         for(let x0=0;x0<w;x0+=step) {
            const values=[];let sx=0,sy=0;
            for(let y=y0;y<Math.min(h,y0+step);++y)
               for(let x=x0;x<Math.min(w,x0+step);++x) {
                  const p=y*w+x;
                  if(!sourceExclusion[p] && (manual ? allowed[p] : !excluded[p]) && isFinite(plane[p])) {
                     values.push(plane[p]);sx+=x;sy+=y;
                  }
               }
            if(values.length>=9) points.push([sx/values.length,sy/values.length,acMedian(values)]);
         }
      }
      return points;
   }
   // Regularized thin plate spline, with bounded spatially distributed donors.
   function model(points,o) {
      const cx=(o.x0+o.x1)/2,cy=(o.y0+o.y1)/2;
      const scale=Math.max(24,(o.x1-o.x0)/2,(o.y1-o.y0)/2);
      const hx=(o.x1-o.x0)/2/scale,hy=(o.y1-o.y0)/2/scale;
      const distance=p=>Math.hypot(Math.max(0,Math.abs(p[0])-hx),Math.max(0,Math.abs(p[1])-hy));
      let samples=points.map(p=>[(p[0]-cx)/scale,(p[1]-cy)/scale,p[2]]);
      samples.sort((a,b)=>distance(a)-distance(b));
      if(!o.regions.length) {
         const radius=Math.max(0.45,distance(samples[Math.min(95,samples.length-1)]));
         samples=samples.filter(p=>distance(p)<=radius);
      }
      // Aggregate nearby tile medians into spatial bins, not just nearest points.
      let cell=0.12, nodes=[];
      do {
         const bins=new Map();
         for(const p of samples) {
            const key=Math.floor(p[0]/cell)+','+Math.floor(p[1]/cell);
            if(!bins.has(key)) bins.set(key,[]);
            bins.get(key).push(p);
         }
         nodes=Array.from(bins.values(),v=>[acMedian(v.map(p=>p[0])),acMedian(v.map(p=>p[1])),acMedian(v.map(p=>p[2]))]);
         cell*=1.4;
      } while(nodes.length>96);
      const n=nodes.length, base=acMedian(nodes.map(p=>p[2]));
      const kernel=r2=>r2>1e-20?0.5*r2*Math.log(r2):0;
      function solve(weights) {
         const m=Array.from({length:n+3},()=>new Float64Array(n+4));
         for(let i=0;i<n;++i) {
            for(let j=0;j<n;++j) {const dx=nodes[i][0]-nodes[j][0],dy=nodes[i][1]-nodes[j][1];m[i][j]=kernel(dx*dx+dy*dy);}
            m[i][i]+=0.035/weights[i];
            m[i][n]=m[n][i]=1;m[i][n+1]=m[n+1][i]=nodes[i][0];m[i][n+2]=m[n+2][i]=nodes[i][1];
            m[i][n+3]=nodes[i][2]-base;
         }
         for(let j=0;j<n+3;++j) {
            let k=j;for(let i=j+1;i<n+3;++i)if(Math.abs(m[i][j])>Math.abs(m[k][j]))k=i;
            [m[j],m[k]]=[m[k],m[j]];const d=m[j][j];if(Math.abs(d)<1e-12)return null;
            for(let l=j;l<n+4;++l)m[j][l]/=d;
            for(let i=0;i<n+3;++i)if(i!==j){const t=m[i][j];for(let l=j;l<n+4;++l)m[i][l]-=t*m[j][l];}
         }
         return m.map(row=>row[n+3]);
      }
      function evalAt(beta,x,y) {
         let v=base+beta[n]+beta[n+1]*x+beta[n+2]*y;
         for(let i=0;i<n;++i){const dx=x-nodes[i][0],dy=y-nodes[i][1];v+=beta[i]*kernel(dx*dx+dy*dy);}
         return v;
      }
      let weights=new Float64Array(n).fill(1),beta=null;
      if(n>=3)for(let iter=0;iter<4;++iter) {
         const next=solve(weights);if(!next)break;beta=next;
         const residual=nodes.map(p=>p[2]-evalAt(beta,p[0],p[1]));
         const med=acMedian(residual.slice()),sigma=Math.max(1e-12,1.4826*acMedian(residual.map(v=>Math.abs(v-med))));
         for(let i=0;i<n;++i)weights[i]=Math.max(0.03,1/Math.max(1,Math.abs(residual[i]-med)/(2.5*sigma)));
      }
      // Evaluate on a small grid and bilinearly interpolate for large masks.
      const step=Math.max(1,Math.ceil(Math.max(o.x1-o.x0,o.y1-o.y0)/160));
      const gw=Math.ceil((o.x1-o.x0)/step)+2,gh=Math.ceil((o.y1-o.y0)/step)+2;
      const grid=new Float64Array(gw*gh);
      for(let j=0;j<gh;++j)for(let i=0;i<gw;++i) {
         const v=beta?evalAt(beta,(o.x0+i*step-cx)/scale,(o.y0+j*step-cy)/scale):base;
         grid[j*gw+i]=Number.isFinite(v)?Math.max(1e-12,v):Math.max(1e-12,base);
      }
      return (x,y)=>{const u=(x-o.x0)/step,v=(y-o.y0)/step,i=Math.floor(u),j=Math.floor(v),a=u-i,b=v-j,k=j*gw+i;return (1-b)*((1-a)*grid[k]+a*grid[k+1])+b*((1-a)*grid[k+gw]+a*grid[k+gw+1]);};
   }
   let automatic=null;
   const trees=objects.map(o=>{
      if(o.regions.length) {
         const allowed=new Uint8Array(mask.length);
         for(const region of o.regions) for(let p=0;p<mask.length;++p) allowed[p]|=region[p];
         const samples=donors(allowed,true);
         if(!samples.length) throw new Error('Background source is too small. Enlarge its lasso.');
         return model(samples,o);
      }
      if(!automatic) {
         const samples=donors(null,false);
         if(!samples.length) throw new Error('Not enough clear background. Mark a source outside the object mask.');
         automatic=samples;
      }
      return model(automatic,o);
   });
   const result=plane.slice();
   // Small spatial blocks limit working memory. There is no copied texture and
   // no feather containing the original galaxy light inside the object mask.
   for(let by=0;by<h;by+=16) {
      tick('Filling object areas — '+Math.round(100*by/h)+'%');
      for(let bx=0;bx<w;bx+=16)
         for(let y=by;y<Math.min(by+16,h);++y)
            for(let x=bx;x<Math.min(bx+16,w);++x) {
               const p=y*w+x;
               if(mask[p]) result[p]=trees[labels[p]](x,y);
            }
   }
   return result;
}

function acMMTResidual(plane,w,h,tick,detailCount=2) {
   if(!Number.isInteger(detailCount) || detailCount<1 || detailCount>6) throw new Error("MMT layers must be 1–6.");
   if(typeof MultiscaleMedianTransform === "undefined")
      throw new Error("PixInsight MultiscaleMedianTransform is not available.");
   const process=new MultiscaleMedianTransform;
   // Reuse the native parameter row schema for this installed PI version.
   const defaults=process.layers;
   if(!defaults || !defaults.length || !Array.isArray(defaults[0]))
      throw new Error("Could not read the native MMT layer parameters.");
   const detail=defaults[0], residual=defaults[defaults.length-1];
   const layers=[];
   for(let i=0;i<=detailCount;++i) {
      const row=(i===detailCount?residual:detail).slice();
      row[0]=i===detailCount; // Only residual R is enabled.
      if(row.length>1) row[1]=false; // No bias.
      if(row.length>2) row[2]=0;
      if(row.length>3) row[3]=false; // No noise reduction.
      layers.push(row);
   }
   process.layers=layers;
   const temporary=acCreateOutput("AstroClean_MMT_Work",w,h,1);
   try {
      acEditOutput(temporary,image=>{
         for(let y=0;y<h;y+=64) {
            const rows=Math.min(64,h-y);
            acWriteRows(image,y,rows,0,plane.subarray(y*w,(y+rows)*w));
         }
      });
      tick("Native MMT — "+detailCount+" detail layers disabled, R only");
      if(!process.executeOn(temporary.mainView,false))
         throw new Error("MultiscaleMedianTransform did not complete.");
      const result=new Float32Array(plane.length);
      for(let y=0;y<h;y+=64) {
         tick();
         const rows=Math.min(64,h-y);
         result.set(acReadRows(temporary.mainView.image,y,rows,0),y*w);
      }
      return result;
   } finally {
      if(!temporary.isNull) temporary.forceClose();
   }
}

// Optional smoothing of the rebuilt flat only, in a narrow band on both sides
// of the object boundary. The binary protection mask is never blurred.
function acSoftenBoundary(rebuilt,mask,w,h,pixels,tick,factor=AC_FLAT_FACTOR) {
   if(!(pixels>0)) return rebuilt;
   const band=Math.max(1,Math.ceil(pixels/factor));
   const dist=new Int32Array(mask.length).fill(-1),queue=new Int32Array(mask.length);
   let tail=0;
   for(let p=0;p<mask.length;++p) {
      let edge=false;acNeighbors(p,w,h,false,q=>{if(mask[q]!==mask[p])edge=true;});
      if(edge){dist[p]=0;queue[tail++]=p;}
   }
   for(let head=0;head<tail;++head) {
      const p=queue[head];if(dist[p]>=band)continue;
      acNeighbors(p,w,h,false,q=>{if(dist[q]<0){dist[q]=dist[p]+1;queue[tail++]=q;}});
   }
   const radius=Math.max(1,Math.round(band/3));
   function box(src,horizontal) {
      const dst=new Float32Array(src.length),outer=horizontal?h:w,length=horizontal?w:h;
      for(let line=0;line<outer;++line) {
         if(line%64===0)tick("Softening flat boundary");
         const at=i=>horizontal?line*w+acReflect(i,w):acReflect(i,h)*w+line;
         let sum=0;for(let k=-radius;k<=radius;++k)sum+=src[at(k)];
         for(let i=0;i<length;++i){dst[at(i)]=sum/(2*radius+1);sum+=src[at(i+radius+1)]-src[at(i-radius)];}
      }
      return dst;
   }
   let smooth=rebuilt;
   for(let k=0;k<3;++k)smooth=box(box(smooth,true),false);
   const out=rebuilt.slice();
   for(let p=0;p<out.length;++p)if(dist[p]>=0 && dist[p]<band){
      const t=dist[p]/band,weight=1-t*t*(3-2*t);
      out[p]=rebuilt[p]*(1-weight)+smooth[p]*weight;
   }
   return out;
}

function acFlat(plane,mask,w,h,tick,regions,detailCount=2,edgePixels=0,factor=AC_FLAT_FACTOR) {
   const rebuilt=acLocalFill(plane,mask,w,h,regions,tick);
   return acMMTResidual(acSoftenBoundary(rebuilt,mask,w,h,edgePixels,tick,factor),w,h,tick,detailCount);
}

// Local additive repair: estimate background on starless data, transfer only
// its smooth difference to the matching image with stars. No flat division.
function acRemoveLocal(original,starless,latest,target,protection,w,h,regions,tick,outputName="AstroClean_Removed") {
   if(!target.some(v=>v)) throw new Error("Draw a lasso around the image feature to remove.");
   const W=original.image.width,H=original.image.height;
   let planes;
   if(latest && latest.repairPlanes) planes=latest.repairPlanes.map(p=>p.slice());
   else {
      planes=acDownsampleImage(starless.image,4,false,tick);
      if(latest && latest.flats) {
         for(let c=0;c<planes.length;++c) {
            const norm=acMedian(acUpsample(latest.flats[c],latest.flatW || w,latest.flatH || h,W,H));
            const smallFlat=acUpsample(latest.flats[c],latest.flatW || w,latest.flatH || h,w,h);
            for(let p=0;p<planes[c].length;++p) planes[c][p]*=norm/smallFlat[p];
         }
      }
   }
   const excluded=Uint8Array.from(target,(v,i)=>v || protection[i] ? 1:0);
   const distance=new Int32Array(target.length),queue=new Int32Array(target.length);
   let tail=0;
   // Distance inside lasso. Border pixels have zero change; full strength 4
   // working pixels (about 16 original pixels) inside the boundary.
   for(let p=0;p<target.length;++p) if(target[p]) {
      let boundary=p%w===0 || p%w===w-1 || p<w || p>=w*(h-1);
      acNeighbors(p,w,h,false,q=>{if(!target[q]) boundary=true;});
      if(boundary) {distance[p]=1;queue[tail++]=p;}
   }
   for(let head=0;head<tail;++head) {
      const p=queue[head];acNeighbors(p,w,h,false,q=>{
         if(target[q] && !distance[q]) {distance[q]=distance[p]+1;queue[tail++]=q;}
      });
   }
   const differences=planes.map((plane,c)=>{
      tick("Local Remove — channel "+(c+1));
      const filled=acAlternativeFill(plane,target,w,h,regions,tick,excluded);
      const delta=new Float32Array(plane.length);
      for(let p=0;p<plane.length;++p) if(target[p]) {
         let t=Math.min(1,Math.max(0,(distance[p]-1)/4));t=t*t*(3-2*t);
         delta[p]=(plane[p]-filled[p])*t;
         plane[p]-=delta[p];
      }
      return delta;
   });
   const output=acCreateOutput(outputName,W,H,planes.length);
   try {
      const input=latest ? latest.corrected.image : original.image;
      acEditOutput(output,image=>{
         for(let c=0;c<planes.length;++c) {
            const full=acUpsample(differences[c],w,h,W,H);
            const support=acUpsample(target,w,h,W,H,true);
            for(let y=0;y<H;y+=64) {
               tick("Applying local Remove");const rows=Math.min(64,H-y);
               const samples=acReadRows(input,y,rows,c);
               for(let i=0;i<samples.length;++i) {
                  const p=y*W+i;if(support[p]) samples[i]-=full[p];
               }
               acWriteRows(image,y,rows,c,samples);
            }
         }
      });
      output.keywords=original.window.keywords;
      if(typeof output.copyAstrometricSolution==="function") output.copyAstrometricSolution(original.window);
      output.mainView.stf=original.stf;
      return {mask:protection,flats:null,w,h,window:output,corrected:output.mainView,
         repairPlanes:planes,repairDifferences:differences,repairTarget:target,isRemoval:true,flatW:latest && (latest.flatW || latest.w),flatH:latest && (latest.flatH || latest.h),baseMask:latest ? (latest.flats ? latest.mask : latest.baseMask) : null,baseW:latest ? (latest.flats ? latest.w : latest.baseW) : null,baseH:latest ? (latest.flats ? latest.h : latest.baseH) : null,baseFlats:latest ? (latest.flats || latest.baseFlats) : null};
   } catch(e) {output.forceClose();throw e;}
}

function acReadRows(image,y,rows,channel) {
   const w = image.width;
   const rect = new Rect(0,y,w,y+rows);
   const buffer = new Float32Array(w*rows);
   try {
      image.getSamples(buffer,rect,channel);
      return buffer;
   } catch (typedArrayError) {
      const samples = new Array(buffer.length);
      image.getSamples(samples,rect,channel);
      return Float32Array.from(samples);
   }
}

function acWriteRows(image,y,rows,channel,buffer) {
   const rect = new Rect(0,y,image.width,y+rows);
   try {
      image.setSamples(buffer,rect,channel);
   } catch (typedArrayError) {
      image.setSamples(Array.from(buffer),rect,channel);
   }
}

function acDownsampleImage(image,factor,luminance,tick) {
   const width = Math.floor(image.width/factor);
   const height = Math.floor(image.height/factor);
   const channels = image.isColor ? 3 : 1;
   const output = Array.from({length:luminance ? 1 : channels},
      () => new Float32Array(width*height));
   const weights = [0.2126,0.7152,0.0722];
   const blockRows = 64;
   for (let start = 0; start < height; start += blockRows/factor) {
      tick();
      const rows = Math.min(height-start,blockRows/factor);
      for (let c = 0; c < channels; ++c) {
         const data = acReadRows(image,start*factor,rows*factor,c);
         const target = output[luminance ? 0 : c];
         const weight = luminance && channels === 3 ? weights[c] : 1;
         for (let y = 0; y < rows; ++y)
            for (let x = 0; x < width; ++x) {
               let sum = 0;
               for (let dy = 0; dy < factor; ++dy) {
                  const source = (y*factor+dy)*image.width+x*factor;
                  for (let dx = 0; dx < factor; ++dx) sum += data[source+dx];
               }
               const i = (start+y)*width+x;
               target[i] += weight*sum/(factor*factor);
               if (!Number.isFinite(target[i]))
                  throw new Error("The starless image contains non-finite samples.");
            }
      }
   }
   return output;
}

function acOutputId(stem) {
   let id = stem;
   for (let next = 2; ImageWindow.windows.some(w => w.mainView.id === id); ++next)
      id = stem + "_" + next;
   return id;
}

function acLoadMask(image,w,h,tick) {
   const xScale = (image.width-1)/(w-1);
   const yScale = (image.height-1)/(h-1);
   const result = new Uint8Array(w*h);
   for (let y = 0; y < h; ++y) {
      if (y % 32 === 0) tick();
      const iy = Math.round(y*yScale);
      const line = acReadRows(image,iy,1,0);
      for (let x = 0; x < w; ++x) {
         const ix = Math.round(x*xScale);
         result[y*w+x] = line[ix] >= 0.5 ? 1 : 0;
      }
   }
   return result;
}

function acCreateOutput(name,W,H,channels) {
   const window = new ImageWindow(W,H,channels,32,true,channels === 3,acOutputId(name));
   if (window.isNull) throw new Error("Could not allocate " + name + ".");
   return window;
}

function acEditOutput(window,write) {
   window.mainView.beginProcess();
   try {
      write(window.mainView.image);
   } finally {
      window.mainView.endProcess();
   }
}

function acShowIntermediates(result,stem) {
   if(!result.flats || result.externalFlat) return;
   const W = result.corrected.image.width, H = result.corrected.image.height;
   const created = [];
   function output(name,count,binary,planeAt) {
      const window = acCreateOutput(name,W,H,count);
      created.push(window);
      acEditOutput(window,image => {
         for (let c = 0; c < count; ++c) {
            const plane = acUpsample(planeAt(c),binary ? result.w : (result.flatW || result.w),binary ? result.h : (result.flatH || result.h),W,H,binary);
            for (let y = 0; y < H; y += 64) {
               const rows = Math.min(64,H-y);
               const block = plane.subarray(y*W,(y+rows)*W);
               acWriteRows(image,y,rows,c,binary ? Float32Array.from(block) : block);
            }
         }
      });
      return window;
   }
   try {
      const mask = output(stem+"_Mask",1,true,() => result.mask);
      const flat = output(stem+"_SyntheticFlat",result.flats.length,false,c => result.flats[c]);
      mask.show();
      flat.show();
   } catch (error) {
      created.reverse().forEach(window => window.forceClose());
      throw error;
   }
}

function acValidateExistingFlat(original,flat,starless=null) {
   if(!original || original.isNull) throw new Error("Select the Original main view.");
   if(!flat || flat.isNull) throw new Error("Select an existing synthetic flat main view.");
   if(original.id===flat.id) throw new Error("Original and synthetic flat must be different views.");
   const a=original.image,b=flat.image;
   if(a.width!==b.width || a.height!==b.height || a.isColor!==b.isColor ||
      a.numberOfChannels!==b.numberOfChannels)
      throw new Error("The flat must have the same dimensions and channels as Original.");
   if(a.width<256 || a.height<256) throw new Error("Both image dimensions must be at least 256 pixels.");
   if(starless && !starless.isNull) {
      if(starless.id===flat.id) throw new Error("Starless must not be the synthetic flat.");
      acValidateInputs(original,starless);
   }
}

function acApplyExistingFlat(original,flat,starless=null) {
   acValidateExistingFlat(original,flat,starless);
   const W=original.image.width,H=original.image.height;
   const channels=original.image.isColor?3:1;
   const progress=new Dialog;progress.windowTitle="AstroClean — Applying existing flat";
   const label=new Label(progress);label.setMinWidth(420);
   progress.sizer=new VerticalSizer;progress.sizer.margin=12;progress.sizer.add(label);
   progress.adjustToContents();progress.show();
   let window=null;
   try {
      const planes=[],norms=[];
      // Read and validate all channels before allocating the result. No smoothing,
      // resampling, mask operation or modification of the supplied flat.
      for(let c=0;c<channels;++c) {
         const plane=new Float32Array(W*H);
         for(let y=0;y<H;y+=64) {
            label.text="Reading flat — channel "+(c+1)+"/"+channels;processEvents();
            const rows=Math.min(64,H-y),block=acReadRows(flat.image,y,rows,c);
            for(const v of block) if(!Number.isFinite(v) || v<=0)
               throw new Error("The flat contains zero, negative or non-finite pixels. Correct these before applying it.");
            plane.set(block,y*W);
         }
         planes.push(plane);norms.push(acMedian(plane));
      }
      window=acCreateOutput("AstroClean_Corrected",W,H,channels);
      acEditOutput(window,image=>{
         for(let c=0;c<channels;++c)
            for(let y=0;y<H;y+=64) {
               label.text="Applying existing flat — channel "+(c+1)+"/"+channels;processEvents();
               const rows=Math.min(64,H-y),samples=acReadRows(original.image,y,rows,c);
               for(let i=0;i<samples.length;++i) samples[i]*=norms[c]/planes[c][y*W+i];
               acWriteRows(image,y,rows,c,samples);
            }
      });
      window.keywords=original.window.keywords;
      if(typeof window.copyAstrometricSolution==="function")window.copyAstrometricSolution(original.window);
      window.mainView.stf=original.stf;
      return {window,corrected:window.mainView,flats:planes,w:W,h:H,
         mask:new Uint8Array(1),externalFlat:true,isRemoval:false};
   } catch(e) {if(window && !window.isNull)window.forceClose();throw e;}
   finally {progress.hide();}
}

function acValidateInputs(original,starless) {
   const input = original.image, background = starless.image;
   const W = input.width, H = input.height;
   if (original.id === starless.id)
      throw new Error("Select two images: the original and its starless counterpart.");
   if (W !== background.width || H !== background.height || input.isColor !== background.isColor)
      throw new Error("The images must have matching dimensions and color spaces.");
   if (W < 256 || H < 256)
      throw new Error("Both image dimensions must be at least 256 pixels.");
}

function acRun(original,starless,editRequest) {
   const factor=editRequest && editRequest.factor || AC_FLAT_FACTOR;
   let aborted=false;
   const progress = new Dialog;
   progress.windowTitle = "AstroClean — Working";
   const progressText = new Label(progress);
   progressText.text = "Preparing background reconstruction…";
   progressText.setMinWidth(480);
   progress.sizer = new VerticalSizer;
   progress.sizer.margin = 14;
   progress.sizer.add(progressText);
   const abortButton=new PushButton(progress);abortButton.text="Abort";
   abortButton.onClick=()=>{aborted=true;abortButton.enabled=false;};
   progress.sizer.add(abortButton);
   progress.adjustToContents();
   const input = original.image, background = starless.image;
   const W = input.width, H = input.height;
   acValidateInputs(original,starless);
   const w = Math.floor(W/4), h = Math.floor(H/4), channels = input.isColor ? 3 : 1;
   const created = [];
   let pulse = 0, stage = "Preparing image", channelStatus = "";
   const tick = text => {
      if (text) stage = text;
      progressText.text = ["◐","◓","◑","◒"][pulse++ % 4] + "  " + channelStatus + stage;
      processEvents();
      if (aborted || console.abortRequested) throw new Error("AstroClean was cancelled.");
   };
   console.abortEnabled = true;
   progress.show();
   const spinner=acSpinner(progressText);
   try {
      console.writeln("<b>AstroClean v2.8.7</b> — object mask");
      let mask;
      if (editRequest && editRequest.preparedMask) {
         mask = editRequest.preparedMask;
      } else if (editRequest) {
         const previous = editRequest.mask.image;
         if (previous.isColor || previous.width !== W || previous.height !== H)
            throw new Error("Select a grayscale AstroClean mask with the same dimensions as the original.");
         mask = acLoadMask(previous,w,h,tick);
      } else {
         const fine = acDownsampleImage(background,2,true,tick)[0];
         mask = acMask(fine,W >> 1,H >> 1,tick);
      }
      const maskCount = mask.reduce((total,v) => total+v,0);
      console.writeln("Masked area: " + (100*maskCount/mask.length).toFixed(2) + "%");
      const flatW=Math.floor(W/factor), flatH=Math.floor(H/factor);
      const fullMask=factor===4 ? mask : acUpsample(mask,w,h,flatW,flatH,true);
      const fullRegions=(editRequest && editRequest.regions || []).map(region=>factor===4 ? region : acUpsample(region,w,h,flatW,flatH,true));
      const planes=acDownsampleImage(background,factor,false,tick);
      const flats = planes.map((plane,c) => {
         channelStatus = "Channel " + (c+1) + "/" + channels + " — ";
         tick("Preparing TPS background");
         const rebuilt=acAlternativeFill(plane,fullMask,flatW,flatH,fullRegions,tick);
         return acMMTResidual(acSoftenBoundary(rebuilt,fullMask,flatW,flatH,editRequest && editRequest.edgePixels || 0,tick,factor),flatW,flatH,tick,editRequest && editRequest.mmtLayers || 2);
      });
      if (flats.some(flat => flat.some(v => !Number.isFinite(v) || v <= 0)))
         throw new Error("The synthetic flat contains nonpositive or non-finite samples.");

      const corrected = acCreateOutput("AstroClean_Corrected",W,H,channels);
      created.push(corrected);
      acEditOutput(corrected,correctedImage => {
         for (let c = 0; c < channels; ++c) {
            channelStatus = "Channel " + (c+1) + "/" + channels + " — ";
            tick("Applying synthetic flat");
            console.writeln("Applying flat, channel " + (c+1));
            const fullFlat = acUpsample(flats[c],flatW,flatH,W,H);
            const normalization = acMedian(fullFlat);
            for (let y = 0; y < H; y += 64) {
               tick();
               const rows = Math.min(64,H-y);
               const values = fullFlat.subarray(y*W,(y+rows)*W);
               const samples = acReadRows(input,y,rows,c);
               for (let i = 0; i < samples.length; ++i)
                  samples[i] *= normalization/values[i];
               acWriteRows(correctedImage,y,rows,c,samples);
            }
         }
      });
      corrected.keywords = original.window.keywords;
      if (typeof corrected.copyAstrometricSolution === "function")
         corrected.copyAstrometricSolution(original.window);
      corrected.mainView.stf = original.stf;
      console.writeln("Done.");
      return {mask,flats,w,h,flatW,flatH,factor,corrected:corrected.mainView,window:corrected,isRemoval:false,edgePixels:editRequest && editRequest.edgePixels || 0,mmtLayers:editRequest && editRequest.mmtLayers || 2,sourceRegions:editRequest && editRequest.regions};
   } catch (error) {
      created.reverse().forEach(window => window.forceClose());
      throw error;
   } finally {
      spinner.stop();progress.hide();
      console.abortEnabled = false;
   }
}

function acLassoMask(mask,w,h,points,canvasWidth,canvasHeight,value) {
   const polygon = points.map(p => ({
      x:p.x*(w-1)/(canvasWidth-1),
      y:p.y*(h-1)/(canvasHeight-1)
   }));
   const revised = mask.slice();
   let changed = 0;
   for (let y = 0; y < h; ++y) {
      const crossings = [];
      for (let i = 0, j = polygon.length-1; i < polygon.length; j = i++) {
         const a = polygon[j], b = polygon[i];
         if ((a.y > y) !== (b.y > y))
            crossings.push(a.x+(y-a.y)*(b.x-a.x)/(b.y-a.y));
      }
      crossings.sort((a,b) => a-b);
      for (let i = 0; i+1 < crossings.length; i += 2) {
         const left = Math.max(0,Math.ceil(crossings[i]));
         const right = Math.min(w,Math.ceil(crossings[i+1]));
         for (let x = left; x < right; ++x) {
            const index = y*w+x;
            if (revised[index] !== value) {
               revised[index] = value;
               ++changed;
            }
         }
      }
   }
   return {mask:revised,changed};
}

function acEditorPixels(view,width,height,factor) {
   const image = view.image;
   const planes = acDownsampleImage(image,factor,false,()=>{});
   const w = Math.floor(image.width/factor), h = Math.floor(image.height/factor);
   return planes.map(plane => acUpsample(plane,w,h,width,height));
}

// Build an exact corrected-starless preview from the same full-resolution gain
// used for the image with stars. Local Remove adds its actual applied delta.
function acCorrectedStarlessPreview(original,starless,result,width,height,factor) {
   if(result.editorStarless) return result.editorStarless;
   const W=original.image.width,H=original.image.height;
   const w=Math.floor(W/factor),h=Math.floor(H/factor);
   const flats=result.flats || result.baseFlats;
   const channels=starless.image.isColor ? 3:1;
   const planes=[];
   for(let c=0;c<channels;++c) {
      const flat=flats ? acUpsample(flats[c],result.flatW || result.w,result.flatH || result.h,W,H) : null;
      const norm=flat ? acMedian(flat):1;
      const small=new Float32Array(w*h);
      for(let y=0;y<h;++y) {
         const top=y*factor;
         const sky=acReadRows(starless.image,top,factor,c);
         const originalRows=result.isRemoval ? acReadRows(original.image,top,factor,c):null;
         const repairedRows=result.isRemoval ? acReadRows(result.corrected.image,top,factor,c):null;
         for(let x=0;x<w;++x) {
            let sum=0;
            for(let dy=0;dy<factor;++dy)
               for(let dx=0;dx<factor;++dx) {
                  const p=dy*W+x*factor+dx;
                  const gain=flat ? norm/flat[top*W+p]:1;
                  let value=sky[p]*gain;
                  if(result.isRemoval) value+=repairedRows[p]-originalRows[p]*gain;
                  sum+=value;
               }
            small[y*w+x]=sum/(factor*factor);
         }
      }
      planes.push(acUpsample(small,w,h,width,height));
   }
   result.editorStarless=planes;
   return planes;
}

function acEditorStretch(planes) {
   return planes.map(plane => {
      const [median,sigma] = acStats(plane);
      const black = median-2.8*sigma;
      return (value,gain,highlight=1) => {
         const scaled = gain*Math.max(0,(value-black)/sigma);
         const knee=3;
         const z=scaled/knee;
         const t=Math.asinh(z);
         const anchor=Math.asinh(1);
         const u=t<=anchor?t:anchor+highlight*(t-anchor);
         return Math.round(255*Math.min(1,u/(u+Math.asinh(25)-anchor)));
      };
   });
}

function acEditorBitmap(planes,stretch,width,height,overlay,gain = 1,highlight=1) {
   const bitmap = new Bitmap(width,height);
   for (let y = 0; y < height; ++y) {
      for (let x = 0; x < width; ++x) {
         const i = y*width+x;
         let r = stretch[0](planes[0][i],gain,highlight);
         let g = r, b = r;
         if (planes.length === 3) {
            g = stretch[1](planes[1][i],gain,highlight);
            b = stretch[2](planes[2][i],gain,highlight);
         }
         if (overlay && overlay[i]) {
            const a=overlay[i];
            r = Math.round((1-0.44*a)*r+20*a);
            g = Math.round((1-0.44*a)*g+112*a);
            b = Math.round((1-0.44*a)*b+48*a);
         }
         bitmap.setPixel(x,y,(0xff000000 | (r<<16) | (g<<8) | b) >>> 0);
      }
   }
   return bitmap;
}

function acMaskEditor(original,starless,first,hidden,resumed) {
   const dialog = new Dialog;
   const external=!!first.externalFlat;
   dialog.windowTitle = external ? "AstroClean — Apply existing flat" : "AstroClean — Mask Editor";
   const source = original.image;
   const scale = Math.min(900/source.width,585/source.height,1);
   const width = Math.max(2,Math.round(source.width*scale));
   const height = Math.max(2,Math.round(source.height*scale));
   // Use a power-of-two sampling factor that divides the block height.
   const factor = Math.min(64,Math.max(1,2**Math.floor(Math.log2(1/scale))));
   const input = acEditorPixels(original,width,height,factor);
   const stretch = acEditorStretch(input);
   let gain = 1, highlight=1;
   let maskVisible=true, appliedMMT=2, appliedEdge=0;
   const canvas = new Control(dialog);
   canvas.setMinSize(400,300);
   canvas.resize(width,height);
   const status = new Label(dialog);
   const instruction = new Label(dialog);
   instruction.text = external ? "Existing flat applied directly. Compare Original and Corrected, then Finish. Flat editing controls are disabled." : "Flat: edit the object mask, then Calculate. Local Remove: draw around the feature, then Remove.";
   instruction.wordWrapping = true;
   status.wordWrapping = true;

   let working = first.mask, applied = first.window ? working : null;
   let history = [], manualEdits=[];
   let sensitivity=50, fineCache=null, appliedFactor=first.factor || AC_FLAT_FACTOR;
   const regions = [];
   let points = [], dragging = false;
   let mode = resumed ? "mask" : "after", busy = false, bitmap;
   let showStarless = !first.window, starlessPixels = null;
   let beforeBitmap = null, afterBitmap = null, flatBitmap = null;
   let latest = first.window ? first : null, previous = null;
   let beforePixels = null, afterPixels = null, flatPlanes = null, flatStretch = null;
   let message = external ? "Existing flat applied — Corrected with stars. Editing controls are locked." : !latest ? "Automatic objects are marked. Add your lasso protection, then Calculate." : resumed ? "Existing mask loaded. Use the lasso to continue editing." :
      "Corrected preview. Choose Mask overlay to edit, or Finish to open the result.";

   function button(label,action) {
      const control = new PushButton(dialog);
      control.text = label;
      control.onClick = action;
      return control;
   }

   function originalPixels() {
      if(!showStarless) return input;
      if(!starlessPixels) starlessPixels=acEditorPixels(starless,width,height,factor);
      return starlessPixels;
   }
   function correctedPixels() {
      if(!latest) return originalPixels();
      if(showStarless) return acCorrectedStarlessPreview(original,starless,latest,width,height,factor);
      if(!afterPixels) afterPixels=acEditorPixels(latest.corrected,width,height,factor);
      return afterPixels;
   }
   function viewStatus() {
      return (mode==="after" ? "Corrected" : "Original") + (showStarless ? " — without stars." : " — with stars.");
   }
   // Visual approximation of the actual two-sided boundary smoothing band.
   function previewMask() {
      const binary=acUpsample(working,first.w,first.h,width,height,true);
      const band=edgeSoft.value*width/source.width;
      if(!band)return binary;
      const d=new Int32Array(binary.length).fill(-1),q=new Int32Array(binary.length);
      let tail=0;
      for(let p=0;p<binary.length;++p) {
         let edge=false;acNeighbors(p,width,height,false,n=>{if(binary[n]!==binary[p])edge=true;});
         if(edge){d[p]=0;q[tail++]=p;}
      }
      for(let head=0;head<tail;++head){const p=q[head];if(d[p]>=band)continue;acNeighbors(p,width,height,false,n=>{if(d[n]<0){d[n]=d[p]+1;q[tail++]=n;}});}
      return Float32Array.from(binary,(v,p)=>{if(d[p]<0)return v;const t=Math.min(1,(d[p]+0.5)/band),a=t*t*(3-2*t);return v?0.5+0.5*a:0.5-0.5*a;});
   }
   function paintBitmap() {
      if (mode === "local") {
         bitmap=acEditorBitmap(correctedPixels(),stretch,width,height,null,gain,highlight);
      } else if (mode === "mask") {
         const overlay = previewMask();
         if (showStarless && !starlessPixels)
            starlessPixels = acEditorPixels(starless,width,height,factor);
         const background = showStarless ? starlessPixels : input;
         bitmap = acEditorBitmap(background,stretch,width,height,maskVisible?overlay:null,gain,highlight);
      } else if (mode === "mask-only") {
         const samples = previewMask();
         bitmap = acEditorBitmap([samples],[value => Math.round(255*value)],width,height);
      } else if (mode === "flat") {
         if (!flatBitmap) {
            if (!flatPlanes) {
               flatPlanes = latest.flats.map(plane =>
                  acUpsample(plane,latest.flatW || latest.w,latest.flatH || latest.h,width,height));
               if(!flatStretch) flatStretch = acEditorStretch(flatPlanes);
            }
            flatBitmap = acEditorBitmap(flatPlanes,flatStretch,width,height,null,gain,highlight);
         }
         bitmap = flatBitmap;
      } else if (mode === "original") {
         bitmap=acEditorBitmap(originalPixels(),stretch,width,height,null,gain,highlight);
      } else {
         bitmap=acEditorBitmap(correctedPixels(),stretch,width,height,null,gain,highlight);
      }
   }

   function refresh(rebuild = false) {
      if (rebuild) paintBitmap();
      add.enabled = remove.enabled = !busy && mode === "mask" && !dragging && points.length >= 3;
      localApply.enabled = !busy && mode === "local" && !dragging && points.length>=3;
      localMode.enabled = !busy && !dragging;
      undoLocal.enabled = !busy && !!latest && !!latest.isRemoval && !!previous;
      visibility.enabled = !busy && mode === "mask";
      visibility.text = maskVisible ? "Mask: On" : "Mask: Off";
      mmt.enabled = !busy && !(latest && latest.isRemoval);
      edgeSoft.enabled=mmt.enabled;
      resolution.enabled=mmt.enabled && !external;
      sensitivityApply.enabled=sensitivitySlider.enabled=!busy && !dragging && !external && !(latest && latest.isRemoval) && !!starless;
      undo.enabled = !busy && history.length > 0;
      recalc.enabled = !busy && !dragging && (working !== applied || mmt.value !== appliedMMT || edgeSoft.value !== appliedEdge || resolutionFactors[resolution.currentItem] !== appliedFactor) && !(latest && latest.isRemoval) &&
         (!!latest || working.some(value => value !== 0));
      toggle.enabled = !busy && !!latest;
      originalView.enabled = !busy && !dragging;
      overlay.enabled = !busy;
      backgroundButton.enabled = !busy && !dragging && !!starless && mode !== "flat" && mode !== "mask-only";
      backgroundButton.text = showStarless ? "Stars: Off" : "Stars: On";
      maskOnly.enabled = !busy && !dragging;
      inspectFlat.enabled = maskOnly.enabled && !!latest && !!latest.flats;
      highlightSlider.enabled = stretchSlider.enabled = !busy && !dragging && mode !== "mask-only";
      finish.enabled = !busy && !!latest && working === applied && (latest && latest.isRemoval || (mmt.value === appliedMMT && edgeSoft.value === appliedEdge && resolutionFactors[resolution.currentItem] === appliedFactor));
      cancel.enabled = !busy;
      toggle.text = "Corrected";
      if(external) for(const control of [overlay,visibility,maskOnly,add,remove,undo,localMode,localApply,undoLocal,recalc,mmt,edgeSoft]) control.enabled=false;
      status.text = message;
      canvas.update();
   }

   const overlay = button("Mask overlay",() => {
      mode = "mask";
      message = "Green pixels are excluded from background estimation.";
      refresh(true);
   });
   const visibility=button("Mask: On",()=>{maskVisible=!maskVisible;refresh(true);});
   const backgroundButton = button("Stars: Off",() => {
      showStarless = !showStarless;
      if(mode==="after" || mode==="original") message=viewStatus();
      refresh(true);
   });
   function applyLasso(value) {
      const polygon=points.slice();
      const change = acLassoMask(working,first.w,first.h,points,width,height,value);
      points = [];
      if (change.changed > 0) {
         history.push({working,manualEdits,sensitivity});
         working = change.mask;
         manualEdits=manualEdits.concat([{polygon,value}]);
         message = "Mask updated. Click Calculate to apply this edit.";
         refresh(true);
      } else {
         message = value ? "The lasso did not cover any unmasked pixels." :
            "The lasso did not cover any masked pixels.";
         refresh();
      }
   }
   const add = button("Add to mask",() => applyLasso(1));
   const remove = button("Remove from mask",() => applyLasso(0));
   const undo = button("Undo mask edit",() => {
      if (!history.length) return;
      const restored = history.pop();
      working = restored.working;
      manualEdits=restored.manualEdits; sensitivity=restored.sensitivity;
      sensitivitySlider.value=sensitivity;sensitivityValue.text=String(sensitivity);
      points = [];
      mode = "mask";
      message = working === applied ? "Showing the last applied mask." : "Mask edit undone.";
      refresh(true);
   });
   const recalc = button("Calculate",() => {
      busy = true;
      message = "Calculating the synthetic flat…";
      refresh();
      try {
         const result = acRun(original,starless,{preparedMask:working,regions,mmtLayers:mmt.value,edgePixels:edgeSoft.value,factor:resolutionFactors[resolution.currentItem]});
         hidden.add(result.window);
         if (previous) {
            previous.window.forceClose();
            hidden.delete(previous.window);
         }
         previous = latest;
         latest = result;
         applied = working;
         appliedMMT=mmt.value;appliedEdge=edgeSoft.value;appliedFactor=resolutionFactors[resolution.currentItem];
         beforePixels = afterPixels;
         afterPixels = flatPlanes = flatStretch = null;
         beforeBitmap = afterBitmap = flatBitmap = null;
         acRememberPair(original,starless);
         points = [];
         mode = "after";
         showStarless=false;
         message = viewStatus();
      } catch (error) {
         message = error.message;
      } finally {
         busy = false;
         refresh(true);
      }
   });
   const localMode=button("Local Remove",()=>{
      points=[];mode="local";
      message="Draw around the feature to remove from the image, then click Remove. This does not edit the object mask.";
      refresh(true);
   });
   const localApply=button("Remove",()=>{
      const target=acLassoMask(new Uint8Array(working.length),first.w,first.h,points,width,height,1).mask;
      busy=true;message="Preparing local Remove…";refresh();processEvents();
      const progress=new Dialog;progress.windowTitle="AstroClean — Local Remove";
      const label=new Label(progress);label.setMinWidth(400);
      progress.sizer=new VerticalSizer;progress.sizer.margin=12;progress.sizer.add(label);
      progress.adjustToContents();progress.show();
      try {
         const tick=text=>{if(text) label.text=text;processEvents();};
         const result=acRemoveLocal(original,starless,latest,target,working,first.w,first.h,regions,tick);
         hidden.add(result.window);
         if(previous){
            previous.window.forceClose();hidden.delete(previous.window);
         }
         if(!latest) {
            const backup=acCreateOutput("AstroClean_BeforeRemove",source.width,source.height,source.isColor?3:1);
            acEditOutput(backup,image=>image.assign(original.image));
            backup.mainView.stf=original.stf;hidden.add(backup);
            latest={window:backup,corrected:backup.mainView,mask:working,w:first.w,h:first.h,flats:null,isRemoval:false};
         }
         previous=latest;latest=result;
         applied=working;appliedMMT=mmt.value;appliedEdge=edgeSoft.value;appliedFactor=resolutionFactors[resolution.currentItem];
         beforePixels=afterPixels=flatPlanes=flatStretch=null;
         beforeBitmap=afterBitmap=flatBitmap=null;points=[];mode="after";
         message="Local Remove applied. Undo Remove restores the previous result.";
         acRememberPair(original,starless);
      } catch(e){message=e.message;}
      finally{progress.hide();busy=false;refresh(true);}
   });
   const undoLocal=button("Undo Remove",()=>{
      if(!latest || !latest.isRemoval || !previous) return;
      latest.window.forceClose();hidden.delete(latest.window);
      latest=previous;previous=null;
      applied=latest.flats ? latest.mask : null;
      appliedMMT=latest.mmtLayers || 2;appliedEdge=latest.edgePixels || 0;appliedFactor=latest.factor || AC_FLAT_FACTOR;
      beforePixels=afterPixels=flatPlanes=flatStretch=null;
      beforeBitmap=afterBitmap=flatBitmap=null;mode="after";
      message="Last local Remove undone.";refresh(true);
   });
   const originalView=button("Original",()=>{
      mode="original";points=[];message=viewStatus();refresh(true);
   });
   const toggle = button("Corrected",() => {
      mode="after";points=[];message=viewStatus();refresh(true);
   });
   const maskOnly = button("Mask only",() => {
      mode = "mask-only";
      message = "White pixels are excluded from background estimation.";
      refresh(true);
   });
   const inspectFlat = button("Synthetic flat",() => {
      mode = "flat";
      message = working === applied ? "Synthetic flat from the latest correction." :
         "Synthetic flat from before the pending mask edits. Click Calculate to update.";
      refresh(true);
   });
   const finish = button("Finish",() => {
      if (!busy && latest && working === applied && (latest.isRemoval || (mmt.value === appliedMMT && edgeSoft.value === appliedEdge && resolutionFactors[resolution.currentItem] === appliedFactor))) dialog.ok();
   });
   const cancel = button("Cancel",() => dialog.cancel());
   const stretchLabel = new Label(dialog);
   const stretchSlider = new Slider(dialog);
   stretchSlider.setRange(0,100);
   stretchSlider.value = 50;
   stretchSlider.toolTip = "Adjust preview brightness without changing the image data.";
   stretchLabel.text = "Background  50";
   stretchSlider.onValueUpdated = value => {
      gain = 2**((value-50)/25);
      stretchLabel.text = "Background  " + value;
      beforeBitmap = afterBitmap = flatBitmap = null;
      refresh(true);
   };
   const stretchRow = new HorizontalSizer;
   stretchRow.spacing = 8;
   stretchRow.add(stretchLabel);
   stretchRow.add(stretchSlider,100);

   const highlightLabel=new Label(dialog);highlightLabel.text="Highlights  50";
   const highlightSlider=new Slider(dialog);highlightSlider.setRange(0,100);highlightSlider.value=50;
   highlightSlider.toolTip="Adjust bright preview structures without modifying linear data.";
   highlightSlider.onValueUpdated=value=>{highlight=2**((value-50)/20);highlightLabel.text="Highlights  "+value;beforeBitmap=afterBitmap=flatBitmap=null;refresh(true);};
   const highlightRow=new HorizontalSizer;highlightRow.spacing=8;highlightRow.add(highlightLabel);highlightRow.add(highlightSlider,100);

   const mmt=new SpinBox(dialog);
   mmt.minValue=1;mmt.maxValue=6;mmt.value=2;
   mmt.toolTip="Disabled MMT detail layers; only R is retained. Applied on the quarter-size model. Default: 2.";
   mmt.onValueUpdated=()=>{message="MMT changed. Click Calculate to apply.";refresh();};
   const mmtLabel=new Label(dialog);mmtLabel.text="MMT layers + R:";
   const edgeSoft=new SpinBox(dialog);edgeSoft.minValue=0;edgeSoft.maxValue=128;edgeSoft.value=0;
   edgeSoft.toolTip="Optional softening of the filled flat boundary, in original-image pixels. 0 preserves the current method. Requires Calculate.";
   edgeSoft.onValueUpdated=()=>{message="Edge softening changed. Click Calculate.";refresh(true);};
   const edgeLabel=new Label(dialog);edgeLabel.text="Edge soften (px):";
   const mmtRow=new HorizontalSizer;mmtRow.spacing=7;mmtRow.add(mmtLabel);mmtRow.add(mmt);mmtRow.addSpacing(14);mmtRow.add(edgeLabel);mmtRow.add(edgeSoft);mmtRow.addStretch();

   const resolutionFactors=[4,2,1];
   const resolution=new ComboBox(dialog);
   for(const label of ["25%","50%","100%"]) resolution.addItem(label);
   resolution.currentItem=resolutionFactors.indexOf(appliedFactor);
   resolution.onItemSelected=()=>{message="Flat resolution changed. Click Calculate.";refresh();};
   const resolutionLabel=new Label(dialog);resolutionLabel.text="Flat resolution:";
   mmtRow.add(resolutionLabel);mmtRow.add(resolution);
   const sensitivitySlider=new Slider(dialog);sensitivitySlider.setRange(0,100);sensitivitySlider.value=50;
   const sensitivityLabel=new Label(dialog);sensitivityLabel.text="Mask sensitivity:";
   const sensitivityValue=new Label(dialog);sensitivityValue.text="50";
   const sensitivityRow=new HorizontalSizer;sensitivityRow.spacing=7;sensitivityRow.add(sensitivityLabel);sensitivityRow.add(sensitivitySlider,100);sensitivityRow.add(sensitivityValue);
   sensitivitySlider.onValueUpdated=v=>{sensitivityValue.text=String(v);};
   function updateSensitivity(){
      const v=sensitivitySlider.value;
      if(busy || v===sensitivity || external || !starless)return;
      busy=true;message="Updating automatic mask…";refresh();processEvents();
      const waiting=new Dialog;
      waiting.windowTitle="AstroClean — Updating mask";
      const waitLabel=new Label(waiting);
      waitLabel.text="Updating automatic mask…";waitLabel.setMinWidth(440);
      waiting.sizer=new VerticalSizer;waiting.sizer.margin=14;waiting.sizer.add(waitLabel);
      waiting.adjustToContents();
      let spinner=null;
      try{
         let pulse=0;
         const tick=()=>{waitLabel.text=["◐","◓","◑","◒"][pulse++%4]+"  Updating automatic mask…";processEvents();};
         waiting.show();tick();spinner=acSpinner(waitLabel);
         if(!fineCache)fineCache=acDownsampleImage(starless.image,2,true,tick)[0];
         let next=acMask(fineCache,original.image.width>>1,original.image.height>>1,tick,v,AC_DEBUG_MASK);
         for(const edit of manualEdits)next=acLassoMask(next,first.w,first.h,edit.polygon,width,height,edit.value).mask;
         history.push({working,manualEdits,sensitivity});
         working=next;sensitivity=v;mode="mask";maskVisible=true;
         message="Mask updated. Click Calculate to apply.";
      }catch(e){sensitivitySlider.value=sensitivity;message=e.message;}
      finally{if(spinner)spinner.stop();waiting.hide();busy=false;refresh(true);}
   }
   // Explicit apply also supports keyboard adjustment and avoids recalculating during drag.
   const sensitivityApply=button("Update mask",updateSensitivity);sensitivityRow.add(sensitivityApply);
   sensitivitySlider.onMouseRelease=()=>updateSensitivity();

   let zoom=1, panX=0, panY=0, panning=false, panStart=null;
   function fit() {
      const scale = Math.min(canvas.width/width,canvas.height/height)*zoom;
      const w = Math.max(1,Math.round(width*scale));
      const h = Math.max(1,Math.round(height*scale));
      return {x:Math.floor((canvas.width-w)/2+panX),y:Math.floor((canvas.height-h)/2+panY),w,h};
   }
   function point(x,y) {
      const area = fit();
      return {
         x:Math.min(width-1,Math.max(0,(x-area.x)*(width-1)/Math.max(1,area.w-1))),
         y:Math.min(height-1,Math.max(0,(y-area.y)*(height-1)/Math.max(1,area.h-1)))
      };
   }
   function clampPan(){
      const a=fit();
      const mx=Math.max(0,(a.w-canvas.width)/2),my=Math.max(0,(a.h-canvas.height)/2);
      panX=Math.max(-mx,Math.min(mx,panX));panY=Math.max(-my,Math.min(my,panY));
   }
   function changeZoom(mult,x,y){
      if(busy || dragging || panning)return;
      const p=point(x,y);zoom=Math.max(1,Math.min(8,zoom*mult));
      const a=fit();panX+=x-(a.x+p.x*(a.w-1)/(width-1));panY+=y-(a.y+p.y*(a.h-1)/(height-1));
      clampPan();canvas.update();
   }
   const zoomRow=new HorizontalSizer;zoomRow.spacing=7;
   zoomRow.add(button("Zoom +",()=>changeZoom(1.5,canvas.width/2,canvas.height/2)));
   zoomRow.add(button("Zoom -",()=>changeZoom(1/1.5,canvas.width/2,canvas.height/2)));
   zoomRow.add(button("Fit",()=>{if(!busy&&!dragging){zoom=1;panX=panY=0;canvas.update();}}));zoomRow.addStretch();
   canvas.onMouseWheel=(x,y,delta)=>{if(delta)changeZoom(delta>0?1.25:0.8,x,y);};
   canvas.onResize = () => {clampPan();canvas.update();};
   canvas.onMousePress = (x,y,mouseButton) => {
      if(!busy && !dragging && (mouseButton===2 || mouseButton===4)){panning=true;panStart={x,y,panX,panY};return;}
      if (busy || (mode !== "mask" && mode !== "local") || mouseButton !== 1) return;
      const area = fit();
      if (x < area.x || y < area.y || x >= area.x+area.w || y >= area.y+area.h) return;
      const start = point(x,y);
      points = [start];
      dragging = true;
      refresh();
   };
   canvas.onMouseMove = (x,y) => {
      if(panning){panX=panStart.panX+x-panStart.x;panY=panStart.panY+y-panStart.y;clampPan();canvas.update();return;}
      if (!dragging) return;
      const last = points[points.length-1];
      const p = point(x,y);
      if (Math.hypot(last.x-p.x,last.y-p.y) >= 2/zoom) points.push(p);
      canvas.update();
   };
   canvas.onMouseRelease = (x,y,mouseButton) => {
      if(panning){if((mouseButton===2 || mouseButton===4) || mouseButton===undefined){panning=false;panStart=null;}return;}
      if (!dragging || (mouseButton!==undefined && mouseButton!==1)) return;
      const p=point(x,y), last=points[points.length-1];
      if(Math.hypot(last.x-p.x,last.y-p.y)>=2/zoom) points.push(p);
      if(points.length<3) points=[];
      dragging = false;
      refresh();
   };
   canvas.onPaint = function() {
      if (!bitmap) return;
      const graphics = new Graphics(this);
      try {
         const area = fit();
         graphics.fillRect(0,0,this.width,this.height,new Brush(0xff181818));
         graphics.drawScaledBitmap(new Rect(area.x,area.y,area.x+area.w,area.y+area.h),bitmap);
         const screen = p => ({
            x:area.x+p.x*(area.w-1)/(width-1),
            y:area.y+p.y*(area.h-1)/(height-1)
         });
         if ((mode === "mask" || mode === "local") && points.length > (dragging ? 1 : 2)) {
            graphics.pen = new Pen(0xffffd262,2);
            for (let i = 1; i < points.length; ++i) {
               const a = screen(points[i-1]), b = screen(points[i]);
               graphics.drawLine(a.x,a.y,b.x,b.y);
            }
            if (!dragging) {
               const a = screen(points[0]), b = screen(points[points.length-1]);
               graphics.drawLine(a.x,a.y,b.x,b.y);
            }
         }
      } finally {
         graphics.end();
      }
   };

   const tools = new HorizontalSizer;
   tools.spacing=7;
   for(const control of [overlay,visibility,backgroundButton,maskOnly]) tools.add(control);
   const edits=new HorizontalSizer;edits.spacing=7;
   for(const control of [add,remove,undo]) edits.add(control);
   const actions=new HorizontalSizer;actions.spacing=7;
   for(const control of [localMode,localApply,undoLocal,recalc]) actions.add(control);
   const outputs=new HorizontalSizer;outputs.spacing=7;
   outputs.add(originalView);outputs.add(toggle);outputs.add(inspectFlat);outputs.add(acHelpButton(dialog));
   outputs.add(acCredit(dialog));outputs.addStretch();outputs.add(cancel);outputs.add(finish);
   dialog.sizer = new VerticalSizer;
   dialog.sizer.margin = 10;
   dialog.sizer.spacing = 7;
   dialog.sizer.add(instruction);
   dialog.sizer.add(canvas,100);
   dialog.sizer.add(zoomRow);
   dialog.sizer.add(status);
   dialog.sizer.add(stretchRow);
   dialog.sizer.add(highlightRow);
   dialog.sizer.add(mmtRow);
   dialog.sizer.add(sensitivityRow);
   dialog.sizer.add(tools);
   dialog.sizer.add(edits);
   dialog.sizer.add(actions);
   dialog.sizer.add(outputs);
   dialog.adjustToContents();
   dialog.resize(960,720);
   refresh(true);
   if (dialog.execute() && working === applied) return latest;
   return null;
}

const acLastPairKey = "AstroClean/LastProcessedPair";

function acViewReference(view) {
   return {id:view.id, path:view.window.filePath || ""};
}

function acReadLastPair() {
   try {
      const text = Settings.read(acLastPairKey,DataType.String);
      const pair = text ? JSON.parse(text) : null;
      return pair && pair.original && pair.starless ? pair : null;
   } catch (error) {
      return null;
   }
}

function acFindOpenView(reference) {
   const views = ImageWindow.windows.map(window => window.mainView);
   if (reference.path) {
      const matches = views.filter(view => view.window.filePath === reference.path);
      return matches.length === 1 ? matches[0] : matches.find(view => view.id === reference.id);
   }
   return views.find(view => view.id === reference.id && !view.window.filePath);
}

function acRememberPair(original,starless) {
   try {
      const pair = {original:acViewReference(original),starless:acViewReference(starless)};
      Settings.write(acLastPairKey,DataType.String,JSON.stringify(pair));
   } catch (error) {
      console.writeln("AstroClean could not remember the image selection: " + error.message);
   }
}

const acHelpText = "AstroClean 2.8.7\n\u00a9 Morten Lerager\n\nIntroduction\nAstroClean removes residual dust shadows, diffuse reflections and other background defects from linear astronomical images using a synthetic flat. The object mask protects genuine astronomical structures from being included in the flat. AstroClean complements, rather than replaces, normal calibration and gradient correction.\n\nRequirements\nUse PixInsight 1.9.6 with V8 and native MultiscaleMedianTransform (MMT). Open a linear Original and a matching linear Starless image with identical geometry and channels. Do not permanently stretch either image. For an existing flat, Starless is optional.\n\nQuick start\n  Select Original and Starless and click  Run AstroClean .  Inspect the automatic green object mask. Add missing galaxies or nebulosity; remove any protected dust shadows or reflections.  Click  Calculate  and inspect Corrected and Synthetic flat.  Adjust the mask or parameters and calculate again if necessary.  Click  Finish  to open the output images. Save them in PixInsight if required.  \n\nMask editor\nGreen areas protect real objects. Use Add to mask or Remove from mask and draw a lasso. Undo mask edit reverses the last edit. Mask On/Off changes overlay visibility only; Mask only displays the protection mask. Mouse wheel zooms; right or middle drag pans; Fit resets the view. Original and Corrected select the before/after image; Stars On/Off switches the corresponding preview. Background and Highlights adjust display stretch only, not the linear data.\n\nParameters\n   Parameter  Purpose    Mask sensitivity (0\u2013100, default 50)  Higher values generally protect more faint structure. Click Update mask to apply and replay manual lasso edits.    MMT layers + R (1\u20136, default 2)  Number of suppressed detail layers in the flat model. Higher values produce a smoother model but can lose small dust shadows.    Edge soften (0\u2013128 px, default 0)  Smooths the filled boundary without changing the binary mask. Try 30\u201350 pixels if a boundary is visible.    Flat resolution (25%, 50%, 100%; default 25%)  Higher resolution can retain smaller defects at greater processing cost.    Click Calculate after changing flat parameters. \n\nSynthetic flat calculation\nProtected areas are reconstructed from surrounding background using regularized thin-plate spline (TPS) interpolation. Boundary softening and MMT filtering form the flat model. The corrected image is computed per channel as Original \u00d7 median(flat) / flat. Dust shadows and reflections should remain visible in the flat, while real astronomical objects should not. Inspect both outputs; a synthetic flat is an estimate, not a measured calibration flat.\n\nLocal Remove\nSelect Local Remove, draw a lasso around an unwanted feature, and click Remove. TPS reconstructs the local background on the starless working image, excluding protected areas from background sampling. The additive difference is transferred to the image with stars, with an approximately 16-pixel internal feather. No synthetic flat is generated or changed by this operation. Undo Remove restores the preceding result (one step). Further removals can be applied successively. Calculate is disabled after a removal; restart with the original inputs for a fresh flat if needed. Real structures within the selection can also be removed; missing stars or hidden detail cannot be reconstructed. Use a small margin around the defect.\n\nResume from mask\nOpen Original, Starless and the saved grayscale mask with matching dimensions. Select Existing mask and click Resume from mask. Edit the protection and calculate a new flat.\n\nApply existing flat\nOpen Original and a matching linear synthetic flat, including one edited externally. Select Existing synthetic flat and click Apply existing flat. The supplied flat is used directly at full resolution, without smoothing or regeneration. Dimensions and channels must match and all flat samples must be finite and positive. Starless is optional. Mask editing, Calculate and Local Remove are disabled. Finish opens only the corrected image.\n\nOutput images\nFinish opens a new image named after Original with the suffix _AstroClean. A calculated mask and synthetic flat are also opened with _Mask and _SyntheticFlat suffixes. If Local Remove follows Calculate, the previously calculated mask and flat are retained; they do not represent the additive local repair. Local Remove alone produces only the repaired image. Cancel discards unpublished results. Source images are not overwritten; save outputs explicitly.\n\nTips and limitations\nIf a dust shadow remains, check that it is not protected by the mask and try a higher flat resolution or fewer MMT layers. If galaxy structure appears in the flat, expand the mask to include faint outer regions. For large objects occupying most of the frame, surrounding background may be insufficient for reliable interpolation. Inspect the result at a consistent stretch. Local Remove is intended for localized defects, not recovery of obscured astronomical detail.\n\nCredits\nCopyright \u00a9 Morten Lerager 2026. The Local Remove workflow was inspired by BlemishBlaster by Franklin Marek; no BlemishBlaster code is embedded.";

function acShowHelp() {
   try {
      if (typeof Dialog.browseScriptDocumentation === "function" && Dialog.browseScriptDocumentation("AstroClean")) return;
   } catch (e) {}
   // Fallback for installations without the external documentation.
   const dialog = new Dialog;
   dialog.windowTitle = "AstroClean — Documentation";
   const text = new TextBox(dialog);
   text.readOnly = true;
   text.text = acHelpText;
   text.setMinWidth(630);
   text.setMinHeight(435);
   const close = new PushButton(dialog);
   close.text = "Close";
   close.onClick = () => dialog.ok();
   const buttons = new HorizontalSizer;
   buttons.addStretch();
   buttons.add(close);
   dialog.sizer = new VerticalSizer;
   dialog.sizer.margin = 12;
   dialog.sizer.spacing = 8;
   dialog.sizer.add(text,100);
   dialog.sizer.add(buttons);
   dialog.adjustToContents();
   dialog.execute();
}

function acHelpButton(parent) {
   let icon = null;
   for (const path of [":/process-interface/browse-documentation.png",":/icons/document.png",":/icons/help.png"]) {
      try {
         const bitmap = new Bitmap(path);
         if (bitmap.width > 0 && bitmap.height > 0) {
            icon = bitmap;
            break;
         }
      } catch (error) {
         // Try the next built-in resource, or use a text button.
      }
   }
   const button = icon ? new ToolButton(parent) : new PushButton(parent);
   if (icon) {
      button.icon = icon;
      button.setFixedSize(26,26);
   } else {
      button.text = "Help";
   }
   button.toolTip = "Browse AstroClean documentation";
   button.onClick = acShowHelp;
   return button;
}

function acCredit(parent) {
   const label = new Label(parent);
   label.text = "\u00a9 Morten Lerager  |  AstroClean v2.8.7";
   return label;
}

function acDialog() {
   console.hide();
   const dialog = new Dialog;
   dialog.windowTitle = "AstroClean v2.8.7";
   const note = new Label(dialog);
   note.text = "Build a flat from Original + Starless, or apply an existing flat to Original. Starless is optional when applying a flat.";
   note.wordWrapping=true;

   function row(title) {
      const label = new Label(dialog);
      label.text = title;
      label.setMinWidth(115);
      const list = new ViewList(dialog);
      list.getMainViews();
      const layout = new HorizontalSizer;
      layout.spacing = 8;
      layout.add(label);
      layout.add(list,100);
      return [layout,list];
   }

   const [originalRow,originalList] = row("Original");
   const [starlessRow,starlessList] = row("Starless");
   const last = acReadLastPair();
   if (last) {
      const original = acFindOpenView(last.original);
      const starless = acFindOpenView(last.starless);
      if (original) originalList.currentView = original;
      if (starless && (!original || starless.id !== original.id))
         starlessList.currentView = starless;
   } else {
      const active = ImageWindow.activeWindow;
      if (!active.isNull) originalList.currentView = active.mainView;
   }

   const [maskRow,maskList] = row("Existing mask");
   const [flatRow,flatList] = row("Existing synthetic flat");
   const editNote = new Label(dialog);
   editNote.text = "Select an existing mask to resume editing.";

   let editRequest = null;
   let manual = false;
   let useExistingFlat=false;
   function accept(action) {
      try {
         const original = originalList.currentView;
         const starless = starlessList.currentView;
         if(action==="flat") {
            acValidateExistingFlat(original,flatList.currentView,starless && !starless.isNull ? starless:null);
            useExistingFlat=true;dialog.ok();return;
         }
         if (!original || original.isNull || !starless || starless.isNull || original.id === starless.id)
            throw new Error("Select two different main views.");
         acValidateInputs(original,starless);
         if (action === "resume") {
            const mask = maskList.currentView;
            if (!mask || mask.isNull)
               throw new Error("Select an existing AstroClean mask.");
            if (mask.image.isColor || mask.image.width !== original.image.width ||
                mask.image.height !== original.image.height)
               throw new Error("Select a grayscale AstroClean mask with the same dimensions as the original.");
            editRequest = {mask,resume:true};
         }
         manual = action === "manual";
         dialog.ok();
      } catch (error) {
         new MessageBox(error.message,"AstroClean",StdIcon_Error,StdButton_Ok).execute();
      }
   }

   const manualStart = new PushButton(dialog);
   manualStart.text = "Run AstroClean";
   manualStart.onClick = () => accept("manual");
   const applyFlat = new PushButton(dialog);
   applyFlat.text="Apply existing flat";
   applyFlat.onClick=()=>accept("flat");
   const resume = new PushButton(dialog);
   resume.text = "Resume from mask";
   resume.onClick = () => accept("resume");
   const cancel = new PushButton(dialog);
   cancel.text = "Cancel";
   cancel.onClick = () => dialog.cancel();
   const buttons = new HorizontalSizer;
   buttons.spacing = 8;
   buttons.add(acHelpButton(dialog));
   buttons.add(acCredit(dialog));
   buttons.addStretch();
   buttons.add(manualStart);
   buttons.add(resume);
   buttons.add(applyFlat);
   buttons.add(cancel);

   dialog.sizer = new VerticalSizer;
   dialog.sizer.margin = 12;
   dialog.sizer.spacing = 8;
   dialog.sizer.add(note);
   dialog.sizer.add(originalRow);
   dialog.sizer.add(starlessRow);
   dialog.sizer.addSpacing(10);
   dialog.sizer.add(editNote);
   dialog.sizer.add(maskRow);
   dialog.sizer.add(flatRow);
   dialog.sizer.add(buttons);
   dialog.adjustToContents();
   dialog.setMinWidth(500);
   if (dialog.execute()) {
      const original = originalList.currentView;
      const starless = starlessList.currentView;
      if(useExistingFlat) {
         const sky=starless && !starless.isNull ? starless:null;
         const appliedFlat=acApplyExistingFlat(original,flatList.currentView,sky);
         const hidden=new Set([appliedFlat.window]);
         try {
            const approved=acMaskEditor(original,sky,appliedFlat,hidden,false);
            if(approved){approved.window.mainView.id=acOutputId(original.id+"_AstroClean");approved.window.show();hidden.delete(approved.window);}
         } finally {for(const window of hidden)if(!window.isNull)window.forceClose();}
         return;
      }
      const w = Math.floor(original.image.width/4), h = Math.floor(original.image.height/4);
      let result;
      if (manual) {
         const waiting = new Dialog;
         waiting.windowTitle = "AstroClean — Finding objects";
         const label = new Label(waiting);
         label.setMinWidth(440);
         waiting.sizer = new VerticalSizer;
         waiting.sizer.margin = 14;
         waiting.sizer.add(label);
         let pulse = 0;
         const tick = () => {
            label.text = ["◐","◓","◑","◒"][pulse++ % 4] + "  Finding objects before lasso editing…";
            processEvents();
         };
         label.text = "Finding objects…";
         waiting.adjustToContents();
         waiting.show();
         try {
            const fine = acDownsampleImage(starless.image,2,true,tick)[0];
            const automatic = acMask(fine,original.image.width >> 1,original.image.height >> 1,tick);
            result = {mask:automatic,w,h};
         } finally {
            waiting.hide();
         }
      } else result = acRun(original,starless,editRequest);
      const hidden = new Set();
      if (result.window) {
         hidden.add(result.window);
         acRememberPair(original,starless);
      }
      let approved = null;
      try {
         approved = acMaskEditor(original,starless,result,hidden,manual || (editRequest && editRequest.resume));
         if (approved) {
            approved.window.mainView.id=acOutputId(original.id+"_AstroClean");
            const intermediates=approved.isRemoval && approved.baseFlats ? Object.assign({},approved,{flats:approved.baseFlats,mask:approved.baseMask,w:approved.baseW,h:approved.baseH}) : approved;
            acShowIntermediates(intermediates,approved.window.mainView.id);
            approved.window.show();
            hidden.delete(approved.window);
         }
      } finally {
         for (const window of hidden)
            if (!window.isNull) window.forceClose();
      }
   }
}

try {
   acDialog();
} catch (error) {
   console.criticalln("AstroClean: " + error.message);
   new MessageBox(error.message,"AstroClean",StdIcon_Error,StdButton_Ok).execute();
}
