#engine v8
#feature-id Lerager > AstroClean
#feature-info AstroClean v2.6.8 experimental. Correct a linear image with a synthetic flat estimated from its starless counterpart.

#include <pjsr/StdButton.jsh>
#include <pjsr/StdIcon.jsh>

"use strict";

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

function acMask(fine,fw,fh,tick) {
   const w = fw >> 1, h = fh >> 1;
   const coarse = acHalf(fine,fw,fh);
   const levels = acPyramid(coarse,w,h,6,tick);
   const extended = levels[2].map((v,i) => v-levels[5][i]);
   const broad = levels[3].map((v,i) => v-levels[6][i]);
   const [em,es] = acStats(extended), [bm,bs] = acStats(broad);
   const support = Uint8Array.from(extended,(v,i) => v > em+1.75*es || broad[i] > bm+2.8*bs);
   const seeds = new Uint8Array(coarse.length);
   acComponents(Uint8Array.from(extended,v => v > em+10*es),w,h,component => {
      if (component.length >= 50 && component.some(i => extended[i] > 18*es))
         for (const i of component) seeds[i] = 1;
   });
   const grown = acFlood(seeds,support,w,h);
   const open = grown.map(v => 1-v);
   const boundary = open.map((v,i) => v && (i < w || i >= w*(h-1) || i%w === 0 || i%w === w-1));
   const closed = acFlood(boundary,open,w,h).map(v => 1-v);
   const large = acDilate(closed,w,h,18);

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
   return acDilate(small,w,h,7).map((v,i) => v || large[i]);
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
function acSoftenBoundary(rebuilt,mask,w,h,pixels,tick) {
   if(!(pixels>0)) return rebuilt;
   const band=Math.max(1,Math.ceil(pixels/4));
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

function acFlat(plane,mask,w,h,tick,regions,detailCount=2,edgePixels=0) {
   const rebuilt=acLocalFill(plane,mask,w,h,regions,tick);
   return acMMTResidual(acSoftenBoundary(rebuilt,mask,w,h,edgePixels,tick),w,h,tick,detailCount);
}

// Local additive repair: estimate background on starless data, transfer only
// its smooth difference to the matching image with stars. No flat division.
function acRemoveLocal(original,starless,latest,target,protection,w,h,regions,tick) {
   if(!target.some(v=>v)) throw new Error("Draw a lasso around the image feature to remove.");
   const W=original.image.width,H=original.image.height;
   let planes;
   if(latest && latest.repairPlanes) planes=latest.repairPlanes.map(p=>p.slice());
   else {
      planes=acDownsampleImage(starless.image,4,false,tick);
      if(latest && latest.flats) {
         for(let c=0;c<planes.length;++c) {
            const norm=acMedian(acUpsample(latest.flats[c],w,h,W,H));
            for(let p=0;p<planes[c].length;++p) planes[c][p]*=norm/latest.flats[c][p];
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
      const filled=acLocalFill(plane,target,w,h,regions,tick,excluded);
      const delta=new Float32Array(plane.length);
      for(let p=0;p<plane.length;++p) if(target[p]) {
         let t=Math.min(1,Math.max(0,(distance[p]-1)/4));t=t*t*(3-2*t);
         delta[p]=(plane[p]-filled[p])*t;
         plane[p]-=delta[p];
      }
      return delta;
   });
   const output=acCreateOutput("AstroClean_Removed",W,H,planes.length);
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
         repairPlanes:planes,isRemoval:true,baseFlats:latest ? (latest.flats || latest.baseFlats) : null};
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

function acShowIntermediates(result) {
   if(!result.flats || result.externalFlat) return;
   const W = result.corrected.image.width, H = result.corrected.image.height;
   const created = [];
   function output(name,count,binary,planeAt) {
      const window = acCreateOutput(name,W,H,count);
      created.push(window);
      acEditOutput(window,image => {
         for (let c = 0; c < count; ++c) {
            const plane = acUpsample(planeAt(c),result.w,result.h,W,H,binary);
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
      const mask = output("AstroClean_Mask",1,true,() => result.mask);
      const flat = output("AstroClean_SyntheticFlat",result.flats.length,false,c => result.flats[c]);
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
   const progress = new Dialog;
   progress.windowTitle = "AstroClean — Working";
   const progressText = new Label(progress);
   progressText.text = "Preparing background reconstruction…";
   progressText.setMinWidth(480);
   progress.sizer = new VerticalSizer;
   progress.sizer.margin = 14;
   progress.sizer.add(progressText);
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
      if (console.abortRequested) throw new Error("AstroClean was cancelled.");
   };
   console.abortEnabled = true;
   progress.show();
   try {
      console.writeln("<b>AstroClean v2.6.8 experimental</b> — object mask");
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
      const planes = acDownsampleImage(background,4,false,tick);
      const flats = planes.map((plane,c) => {
         channelStatus = "Channel " + (c+1) + "/" + channels + " — ";
         tick("Preparing background");
         console.writeln("Reconstructing background, channel " + (c+1));
         return acFlat(plane,mask,w,h,tick,editRequest && editRequest.regions,editRequest && editRequest.mmtLayers || 2,editRequest && editRequest.edgePixels || 0);
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
            const fullFlat = acUpsample(flats[c],w,h,W,H);
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
      return {mask,flats,w,h,corrected:corrected.mainView,window:corrected,isRemoval:false,edgePixels:editRequest && editRequest.edgePixels || 0,mmtLayers:editRequest && editRequest.mmtLayers || 2,sourceRegions:editRequest && editRequest.regions};
   } catch (error) {
      created.reverse().forEach(window => window.forceClose());
      throw error;
   } finally {
      progress.hide();
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
      const flat=flats ? acUpsample(flats[c],result.w,result.h,W,H) : null;
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
      return (value,gain) => {
         const scaled = gain*Math.max(0,(value-black)/sigma);
         return Math.round(255*Math.min(1,Math.asinh(scaled/3)/Math.asinh(25)));
      };
   });
}

function acEditorBitmap(planes,stretch,width,height,overlay,gain = 1) {
   const bitmap = new Bitmap(width,height);
   for (let y = 0; y < height; ++y) {
      for (let x = 0; x < width; ++x) {
         const i = y*width+x;
         let r = stretch[0](planes[0][i],gain);
         let g = r, b = r;
         if (planes.length === 3) {
            g = stretch[1](planes[1][i],gain);
            b = stretch[2](planes[2][i],gain);
         }
         if (overlay && overlay[i]) {
            r = Math.round(0.56*r+20);
            g = Math.round(0.56*g+112);
            b = Math.round(0.56*b+48);
         }
         bitmap.setPixel(x,y,(0xff000000 | (r<<16) | (g<<8) | b) >>> 0);
      }
   }
   return bitmap;
}

// Inward feather on the flat grid; source planes remain immutable.
function acAdjustFlatPlanes(planes,weight,amount) {
   const factor=Math.pow(2,amount/50);
   return planes.map(p=>Float32Array.from(p,(v,i)=>v*(1+weight[i]*(factor-1))));
}
function acFlatSelectionWeight(mask,w,h,pixels,scale) {
   const d=new Int32Array(mask.length),q=[];
   for(let y=0;y<h;++y) for(let x=0;x<w;++x) {
      const i=y*w+x;
      if(mask[i] && (!x || !y || x===w-1 || y===h-1 || !mask[i-1] || !mask[i+1] || !mask[i-w] || !mask[i+w])) {d[i]=1;q.push(i);}
   }
   for(let k=0;k<q.length;++k) {
      const i=q[k],x=i%w,y=Math.floor(i/w);
      for(const j of [x?i-1:-1,x<w-1?i+1:-1,y?i-w:-1,y<h-1?i+w:-1])
         if(j>=0 && mask[j] && !d[j]) {d[j]=d[i]+1;q.push(j);}
   }
   return Float32Array.from(mask,(v,i)=>{const t=pixels>0?Math.min(1,Math.max(0,(d[i]-1)*scale/pixels)):1;return v*t*t*(3-2*t);});
}
function acApplyAdjustedFlat(original,base,planes) {
   const W=original.image.width,H=original.image.height;
   let window=null;
   try {
      window=acCreateOutput("AstroClean_Corrected",W,H,planes.length);
      acEditOutput(window,image=>{
         for(let c=0;c<planes.length;++c) {
            const full=acUpsample(planes[c],base.w,base.h,W,H),norm=acMedian(full);
            for(let y=0;y<H;y+=64) {
               const rows=Math.min(64,H-y),samples=acReadRows(original.image,y,rows,c);
               for(let i=0;i<samples.length;++i) samples[i]*=norm/full[y*W+i];
               acWriteRows(image,y,rows,c,samples);processEvents();
            }
         }
      });
      window.keywords=original.window.keywords;
      if(typeof window.copyAstrometricSolution==="function")window.copyAstrometricSolution(original.window);
      window.mainView.stf=original.stf;
      return Object.assign({},base,{window,corrected:window.mainView,flats:planes,editorStarless:null});
   } catch(e) {if(window)window.forceClose();throw e;}
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
   let gain = 1;
   let flatEditing=false,flatWeight=null,flatDraft=null,flatUndo=null;
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
   let history = [];
   let regions = [], appliedRegions = regions;
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
      if(flatDraft) {
         const base=showStarless ? acCorrectedStarlessPreview(original,starless,latest,width,height,factor) : (afterPixels || (afterPixels=acEditorPixels(latest.corrected,width,height,factor)));
         return base.map((p,c)=>{
            const a=acUpsample(latest.flats[c],latest.w,latest.h,width,height);
            const b=acUpsample(flatDraft[c],latest.w,latest.h,width,height);
            const ratio=acMedian(flatDraft[c])/acMedian(latest.flats[c]);
            return Float32Array.from(p,(v,i)=>v*ratio*a[i]/b[i]);
         });
      }
      if(showStarless) return acCorrectedStarlessPreview(original,starless,latest,width,height,factor);
      if(!afterPixels) afterPixels=acEditorPixels(latest.corrected,width,height,factor);
      return afterPixels;
   }
   function viewStatus() {
      return (mode==="after" ? "Corrected" : "Original") + (showStarless ? " — without stars." : " — with stars.");
   }
   function paintBitmap() {
      if (mode === "local") {
         bitmap=acEditorBitmap(correctedPixels(),stretch,width,height,null,gain);
         const sources=new Uint8Array(working.length);
         for(const region of regions) for(let i=0;i<region.length;++i) sources[i]|=region[i];
         const sky=acUpsample(sources,first.w,first.h,width,height,true);
         for(let y=0;y<height;++y) for(let x=0;x<width;++x)
            if(sky[y*width+x]) bitmap.setPixel(x,y,0xff5599dd);
      } else if (mode === "mask") {
         const overlay = acUpsample(working,first.w,first.h,width,height,true);
         if (showStarless && !starlessPixels)
            starlessPixels = acEditorPixels(starless,width,height,factor);
         const background = showStarless ? starlessPixels : input;
         const skyOverlay = new Uint8Array(working.length);
         for (const region of regions)
            for (let i = 0; i < region.length; ++i) skyOverlay[i] |= region[i];
         const sky = acUpsample(skyOverlay,first.w,first.h,width,height,true);
         bitmap = acEditorBitmap(background,stretch,width,height,maskVisible?overlay:null,gain);
         for (let y = 0; y < height; ++y)
            for (let x = 0; x < width; ++x)
               if (sky[y*width+x]) bitmap.setPixel(x,y,0xff5599dd);
      } else if (mode === "mask-only") {
         const samples = acUpsample(working,first.w,first.h,width,height,true);
         bitmap = acEditorBitmap([samples],[value => value ? 255 : 0],width,height);
      } else if (mode === "flat") {
         if (!flatBitmap) {
            if (!flatPlanes) {
               flatPlanes = (flatDraft || latest.flats).map(plane =>
                  acUpsample(plane,latest.w,latest.h,width,height));
               if(!flatStretch) flatStretch = acEditorStretch(flatPlanes);
            }
            flatBitmap = acEditorBitmap(flatPlanes,flatStretch,width,height,null,gain);
         }
         bitmap = flatBitmap;
      } else if (mode === "original") {
         bitmap=acEditorBitmap(originalPixels(),stretch,width,height,null,gain);
      } else {
         bitmap=acEditorBitmap(correctedPixels(),stretch,width,height,null,gain);
      }
   }

   function refresh(rebuild = false) {
      if (rebuild) paintBitmap();
      add.enabled = remove.enabled = !busy && mode === "mask" && !dragging && points.length >= 3;
      markSky.enabled = !busy && (mode === "mask" || mode === "local") && !dragging && points.length>=3;
      localApply.enabled = !busy && mode === "local" && !dragging && points.length>=3;
      localMode.enabled = !busy && !dragging;
      undoLocal.enabled = !busy && !!latest && !!latest.isRemoval && !!previous;
      visibility.enabled = !busy && mode === "mask";
      visibility.text = maskVisible ? "Mask: On" : "Mask: Off";
      mmt.enabled = !busy && !(latest && latest.isRemoval);
      edgeSoft.enabled=mmt.enabled;
      undo.enabled = !busy && history.length > 0;
      recalc.enabled = !busy && !dragging && (working !== applied || regions !== appliedRegions || mmt.value !== appliedMMT || edgeSoft.value !== appliedEdge) && !(latest && latest.isRemoval) &&
         (!!latest || working.some(value => value !== 0));
      toggle.enabled = !busy && !!latest;
      originalView.enabled = !busy && !dragging;
      overlay.enabled = !busy;
      backgroundButton.enabled = !busy && !dragging && !!starless && mode !== "flat" && mode !== "mask-only";
      backgroundButton.text = showStarless ? "Stars: Off" : "Stars: On";
      maskOnly.enabled = !busy && !dragging;
      inspectFlat.enabled = maskOnly.enabled && !!latest && !!latest.flats;
      stretchSlider.enabled = !busy && !dragging && mode !== "mask-only";
      finish.enabled = !busy && !!latest && working === applied && regions === appliedRegions && (latest && latest.isRemoval || (mmt.value === appliedMMT && edgeSoft.value === appliedEdge));
      cancel.enabled = !busy;
      toggle.text = "Corrected";
      if(external) for(const control of [overlay,visibility,maskOnly,add,remove,markSky,undo,localMode,localApply,undoLocal,recalc,mmt,edgeSoft]) control.enabled=false;
      flatEdit.enabled=!busy && !external && !!latest && !!latest.flats && working===applied && regions===appliedRegions && mmt.value===appliedMMT && edgeSoft.value===appliedEdge;
      flatSelect.enabled=flatEdit.enabled && flatEditing && !dragging && points.length>=3;
      flatAmount.enabled=flatEdit.enabled && !!flatWeight;
      flatFeather.enabled=flatEdit.enabled && flatEditing && !flatWeight;
      flatApply.enabled=flatEdit.enabled && !!flatDraft;
      flatRevert.enabled=!busy && (!!flatUndo || flatEditing);
      if(flatEditing) {
         for(const c of [overlay,visibility,maskOnly,add,remove,markSky,undo,localMode,localApply,undoLocal,recalc,mmt,edgeSoft,originalView,finish]) c.enabled=false;
      }
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
      const change = acLassoMask(working,first.w,first.h,points,width,height,value);
      points = [];
      if (change.changed > 0) {
         history.push({working,regions});
         working = change.mask;
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
   const markSky = button("Mark background",() => {
      const region = acLassoMask(new Uint8Array(working.length),first.w,first.h,points,width,height,1).mask;
      let count = 0;
      for (let i = 0; i < region.length; ++i)
         if (region[i] && working[i]) {
            message = "Sky selection overlaps protected pixels. Choose a clear area outside the green mask.";
            refresh();
            return;
         } else count += region[i];
      if (count < 9) {
         message = "Choose a larger sky region.";
         refresh();
         return;
      }
      history.push({working,regions});
      regions = regions.concat([region]);
      points = [];
      message = regions.length + " source regions selected. Each source is assigned to the nearest masked object.";
      refresh(true);
   });
   const undo = button("Undo mask edit",() => {
      if (!history.length) return;
      const restored = history.pop();
      working = restored.working;
      regions = restored.regions;
      points = [];
      mode = "mask";
      message = working === applied ? "Showing the last applied mask." : "Mask edit undone.";
      refresh(true);
   });
   const recalc = button("Calculate",() => {
      if(flatUndo){hidden.delete(flatUndo.window);flatUndo.window.forceClose();flatUndo=null;}
      busy = true;
      message = "Calculating the synthetic flat…";
      refresh();
      try {
         const result = acRun(original,starless,{preparedMask:working,regions,mmtLayers:mmt.value,edgePixels:edgeSoft.value});
         hidden.add(result.window);
         if (previous) {
            previous.window.forceClose();
            hidden.delete(previous.window);
         }
         previous = latest;
         latest = result;
         applied = working;
         appliedRegions = regions;
         appliedMMT=mmt.value;appliedEdge=edgeSoft.value;
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
      if(flatUndo){hidden.delete(flatUndo.window);flatUndo.window.forceClose();flatUndo=null;}
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
         if(previous){previous.window.forceClose();hidden.delete(previous.window);}
         if(!latest) {
            const backup=acCreateOutput("AstroClean_BeforeRemove",source.width,source.height,source.isColor?3:1);
            acEditOutput(backup,image=>image.assign(original.image));
            backup.mainView.stf=original.stf;hidden.add(backup);
            latest={window:backup,corrected:backup.mainView,mask:working,w:first.w,h:first.h,flats:null,isRemoval:false};
         }
         previous=latest;latest=result;
         applied=working;appliedRegions=regions;appliedMMT=mmt.value;appliedEdge=edgeSoft.value;
         beforePixels=afterPixels=flatPlanes=flatStretch=null;
         beforeBitmap=afterBitmap=flatBitmap=null;points=[];mode="after";
         message="Local Remove applied to the image. Undo Remove restores the previous result. No synthetic flat was applied by Remove.";
         acRememberPair(original,starless);
      } catch(e){message=e.message;}
      finally{progress.hide();busy=false;refresh(true);}
   });
   const undoLocal=button("Undo Remove",()=>{
      if(!latest || !latest.isRemoval || !previous) return;
      latest.window.forceClose();hidden.delete(latest.window);
      latest=previous;previous=null;
      applied=latest.flats ? latest.mask : null;
      appliedRegions=latest.sourceRegions || regions;
      appliedMMT=latest.mmtLayers || 2;appliedEdge=latest.edgePixels || 0;
      beforePixels=afterPixels=flatPlanes=flatStretch=null;
      beforeBitmap=afterBitmap=flatBitmap=null;mode="after";
      message="Last local Remove undone.";refresh(true);
   });
   const originalView=button("Original",()=>{
      mode="original";points=[];message=viewStatus();refresh(true);
   });
   const toggle = button("Corrected",() => {
      mode="after";if(!flatEditing) points=[];message=viewStatus();refresh(true);
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
      if (!busy && latest && working === applied && regions === appliedRegions && (latest.isRemoval || (mmt.value === appliedMMT && edgeSoft.value === appliedEdge))) dialog.ok();
   });
   const cancel = button("Cancel",() => dialog.cancel());
   const stretchLabel = new Label(dialog);
   const stretchSlider = new Slider(dialog);
   stretchSlider.setRange(-20,20);
   stretchSlider.value = 0;
   stretchSlider.toolTip = "Adjust preview brightness without changing the image data.";
   stretchLabel.text = "Preview stretch  1.00×";
   stretchSlider.onValueUpdated = value => {
      gain = 2**(value/10);
      stretchLabel.text = "Preview stretch  " + gain.toFixed(2) + "×";
      beforeBitmap = afterBitmap = flatBitmap = null;
      refresh(true);
   };
   const stretchRow = new HorizontalSizer;
   stretchRow.spacing = 8;
   stretchRow.add(stretchLabel);
   stretchRow.add(stretchSlider,100);

   const mmt=new SpinBox(dialog);
   mmt.minValue=1;mmt.maxValue=6;mmt.value=2;
   mmt.toolTip="Disabled MMT detail layers; only R is retained. Applied on the quarter-size model. Default: 2.";
   mmt.onValueUpdated=()=>{message="MMT changed. Click Calculate to apply.";refresh();};
   const mmtLabel=new Label(dialog);mmtLabel.text="MMT layers + R:";
   const edgeSoft=new SpinBox(dialog);edgeSoft.minValue=0;edgeSoft.maxValue=128;edgeSoft.value=0;
   edgeSoft.toolTip="Optional softening of the filled flat boundary, in original-image pixels. 0 preserves the current method. Requires Calculate.";
   edgeSoft.onValueUpdated=()=>{message="Edge softening changed. Click Calculate.";refresh();};
   const edgeLabel=new Label(dialog);edgeLabel.text="Edge soften (px):";
   const mmtRow=new HorizontalSizer;mmtRow.spacing=7;mmtRow.add(mmtLabel);mmtRow.add(mmt);mmtRow.addSpacing(14);mmtRow.add(edgeLabel);mmtRow.add(edgeSoft);mmtRow.addStretch();

   function fit() {
      const scale = Math.min(canvas.width/width,canvas.height/height);
      const w = Math.max(1,Math.round(width*scale));
      const h = Math.max(1,Math.round(height*scale));
      return {x:Math.floor((canvas.width-w)/2),y:Math.floor((canvas.height-h)/2),w,h};
   }
   function point(x,y) {
      const area = fit();
      return {
         x:Math.min(width-1,Math.max(0,(x-area.x)*(width-1)/Math.max(1,area.w-1))),
         y:Math.min(height-1,Math.max(0,(y-area.y)*(height-1)/Math.max(1,area.h-1)))
      };
   }
   canvas.onResize = () => canvas.update();
   canvas.onMousePress = (x,y,mouseButton) => {
      if (busy || (mode !== "mask" && mode !== "local" && !(flatEditing && (mode==="after" || mode==="flat"))) || mouseButton !== 1) return;
      const area = fit();
      if (x < area.x || y < area.y || x >= area.x+area.w || y >= area.y+area.h) return;
      points = [point(x,y)];
      dragging = true;
      refresh();
   };
   canvas.onMouseMove = (x,y) => {
      if (!dragging) return;
      const last = points[points.length-1];
      const p = point(x,y);
      if (Math.hypot(last.x-p.x,last.y-p.y) >= 2) points.push(p);
      canvas.update();
   };
   canvas.onMouseRelease = (x,y) => {
      if (!dragging) return;
      points.push(point(x,y));
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
         if ((mode === "mask" || mode === "local" || flatEditing) && points.length > 1) {
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

   function clearFlatPreview() {
      flatDraft=null;flatWeight=null;flatAmount.value=1000;flatValue.text="Flat brightness: 1.0000×";points=[];
      flatPlanes=flatBitmap=null;
   }
   const flatEdit=button("Edit flat",()=>{
      clearFlatPreview();flatEditing=true;
      if(mode!=="flat")mode="after";
      message="Draw on Corrected or Synthetic flat, then Select area. Flat brightness changes both views.";refresh(true);
   });
   const flatSelect=button("Select area",()=>{
      const mask=acLassoMask(new Uint8Array(latest.w*latest.h),latest.w,latest.h,points,width,height,1).mask;
      if(!mask.some(v=>v))return;
      flatWeight=acFlatSelectionWeight(mask,latest.w,latest.h,flatFeather.value,source.width/latest.w);
      message="Adjust Flat brightness: darker flat brightens Corrected. Apply edit commits the result.";refresh();
   });
   const flatAmount=new Slider(dialog);flatAmount.setRange(0,2000);flatAmount.value=1000;
   flatAmount.toolTip="Local flat brightness: 0.9000–1.1000×; centre = 1.0000×. Arrow keys adjust by 0.0001. Preview is approximate; Apply edit computes full resolution.";
   const flatValue=new Label(dialog);flatValue.text="Flat brightness: 1.0000×";
   flatAmount.onValueUpdated=position=>{
      const multiplier=0.9+position/10000;
      const value=50*Math.log2(multiplier);
      if(!flatWeight)return;
      flatDraft=acAdjustFlatPlanes(latest.flats,flatWeight,value);
      flatValue.text="Flat brightness: "+multiplier.toFixed(4)+"×";
      flatPlanes=flatBitmap=null;refresh(true);
   };
   const flatFeather=new SpinBox(dialog);flatFeather.minValue=0;flatFeather.maxValue=512;flatFeather.value=30;
   flatFeather.toolTip="Inward feather in original-image pixels; choose before Select area.";
   const flatFeatherLabel=new Label(dialog);flatFeatherLabel.text="Edit feather (px):";
   const flatApply=button("Apply edit",()=>{
      if(!flatDraft)return;
      busy=true;message="Applying flat edit at full resolution…";refresh();processEvents();
      try {
         const result=acApplyAdjustedFlat(original,latest,flatDraft);hidden.add(result.window);
         if(flatUndo){flatUndo.window.forceClose();hidden.delete(flatUndo.window);}
         flatUndo=latest;latest=result;
         flatEditing=false;clearFlatPreview();afterPixels=null;flatStretch=null;
         mode="after";message="Flat edit applied. Corrected and exported flat use the edited model.";
      } catch(e){message=e.message;}
      finally{busy=false;refresh(true);}
   });
   const flatRevert=button("Undo flat edit",()=>{
      if(flatEditing){flatEditing=false;clearFlatPreview();message="Unapplied flat edit discarded.";}
      else if(flatUndo){latest.window.forceClose();hidden.delete(latest.window);latest=flatUndo;flatUndo=null;clearFlatPreview();afterPixels=null;flatStretch=null;message="Last flat edit undone.";}
      mode="after";refresh(true);
   });
   const flatRow=new HorizontalSizer;flatRow.spacing=7;
   for(const c of [flatEdit,flatSelect,flatFeatherLabel,flatFeather,flatApply,flatRevert])flatRow.add(c);
   const flatSliderRow=new HorizontalSizer;flatSliderRow.spacing=7;flatSliderRow.add(flatValue);flatSliderRow.add(flatAmount,100);

   const tools = new HorizontalSizer;
   tools.spacing=7;
   for(const control of [overlay,visibility,backgroundButton,maskOnly]) tools.add(control);
   const edits=new HorizontalSizer;edits.spacing=7;
   for(const control of [add,remove,markSky,undo]) edits.add(control);
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
   dialog.sizer.add(status);
   dialog.sizer.add(stretchRow);
   dialog.sizer.add(mmtRow);
   dialog.sizer.add(tools);
   dialog.sizer.add(edits);
   dialog.sizer.add(actions);
   dialog.sizer.add(flatRow);
   dialog.sizer.add(flatSliderRow);
   dialog.sizer.add(outputs);
   dialog.adjustToContents();
   dialog.resize(960,720);
   refresh(true);
   if (dialog.execute() && working === applied && regions === appliedRegions) return latest;
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

const acHelpText = "2.6.8: Edit flat — draw on Corrected or Synthetic flat, set Edit feather, Select area, adjust Flat brightness, then Apply edit. Darker flat brightens Corrected. Slider preview is approximate; Apply edit uses full resolution. Undo flat edit discards a pending edit or restores the last applied edit. Finish exports the edited flat. Calculate rebuilds it.\n\nAstroClean 2.6.8 test\n\u00a9 Morten Lerager\n\nRun AstroClean opens automatic object detection and the editor before calculation. Resume from mask reopens saved protection.\n\nSynthetic flat\nGreen marks objects excluded from the flat. Add to mask and Remove from mask edit that protection only. Mask: On/Off changes visibility without changing pixels. MMT layers + R ranges from 1 to 6, default 2, at quarter-size working resolution. Calculate applies the mask, sources and MMT selection. A changed MMT value requires Calculate before Finish.\n\nLocal Remove (experimental)\nClick Local Remove, draw a lasso around the unwanted reflection or other image feature, then click Remove. This action changes the result image directly and does not change the object mask. It works without calculating a synthetic flat first. The source Original and Starless remain untouched. Finish opens a new result.\n\nRemove estimates a local background from starless data, feathers the correction inside the lasso over approximately 16 original pixels, and subtracts only that difference from the matching image with stars. It does not copy background texture onto stars. Faint real structure inside the selection can also be removed. Select only the unwanted feature, with a little margin for feathering.\n\nOptional Mark background selections are clean local sources; avoid protected objects and the feature being removed. Automatic donors exclude the object mask. Undo Remove restores the immediately previous result, one step. Additional Remove operations can be applied successively. Calculate is disabled after Remove; Undo Remove returns to the prior calculated result. For a fresh flat after multiple removes, restart with the untouched inputs. Synthetic flat is disabled after Remove because local repairs are not represented by a flat.\n\nThe editor starts at approximately 4:3, remains resizable and fits the image without stretching its aspect ratio. Tool buttons are grouped on shorter rows.\n\nFinish exports the mask and flat for a flat calculation. For Local Remove it opens only the directly repaired result. Cancel discards unpublished results. Sources and local repair selections are session-only.\n\nInstall AstroClean_v268.js via Execute Script File or Feature Scripts / MyScript. Requires PixInsight 1.9.6 / V8 and native MMT.\n\nValidation: syntax, numerical local-repair tests and mocked editor/process integration. This environment cannot execute native PixInsight; test dialog layout and the experimental Remove on copies. Large-object reconstruction remains outside this version's intended use.\n\nThe local-fill code is independently implemented. BlemishBlaster by Franklin Marek inspired the workflow; no BlemishBlaster code is embedded.\n\n2.6.8 fix: Calculate completion now assigns a strict boolean to Undo Remove availability. Standard flat results explicitly identify themselves as non-removal results. Image-processing calculations are unchanged.\n\n2.6.8 \u2014 Four image views\nOriginal and Corrected select the before/after image. Stars: On/Off switches between with-stars and starless in either view. All four share the same fixed preview stretch and brightness slider. Calculate always ends on Corrected and retains the selected Stars setting.\n\nCorrected Starless is computed from the existing starless input using exactly the same full-resolution flat gain as Corrected. No new star removal is run. For Local Remove, the applied local difference is transferred to the corrected starless preview too. The first corrected-starless display can take extra time; it is cached for that result. Exports still contain the corrected image with stars.\n\nMask overlay remains the input image for editing protection. Mask: Off only hides the overlay; use Corrected to inspect the output.\n\n2.6.8\nCalculate now always ends on Corrected with Stars: On.\nEdge soften (px) is optional, 0\u2013128 original-image pixels, default 0. Start with 32 if the flat has a sharp boundary. Calculate applies the choice. It smooths only a narrow band around the filled flat boundary before MMT. It does not change the binary object mask, and it never directly mixes the original galaxy back into the fill. Nearby unmasked signal can still influence smoothing, so inspect the result. Return to 0 to restore the previous filtering exactly. MMT remains 2 by default.\n\nApply existing flat \u2014 2.6.8\nOpen Original and your edited synthetic flat in PixInsight. Select them in Original and Existing synthetic flat, then click Apply existing flat. Starless is optional in this mode; leave it unselected if you only need views with stars. If selected, it must be a matching starless copy of Original.\n\nThe supplied flat is read at full resolution and used directly with Original \u00d7 median(flat) / flat, independently per channel. It is not smoothed, regenerated or rescaled, and the source views are not modified. Dimensions and channels must match. All flat pixels must be finite and strictly positive.\n\nThe editor opens on Corrected with Stars On. Original/Corrected, Stars On/Off (when Starless was supplied), Synthetic flat and preview stretch remain available. Mask editing, MMT, Edge soften, Calculate and Local Remove are disabled in this mode. An existing mask is ignored.\n\nFinish opens only the new corrected image. It does not duplicate or overwrite the supplied flat. Cancel discards the new result. Save images explicitly in PixInsight if you want files on disk.\n";

function acShowHelp() {
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
   label.text = "\u00a9 Morten Lerager  |  AstroClean v2.6.8";
   return label;
}

function acDialog() {
   console.hide();
   const dialog = new Dialog;
   dialog.windowTitle = "AstroClean v2.6.8 experimental";
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
            if(approved){approved.window.show();hidden.delete(approved.window);}
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
            acShowIntermediates(approved);
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
