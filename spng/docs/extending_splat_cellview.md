### Introduction
This page describes a 'case-study' in extending the a 'complete' graph. The goal here is to lay out my own thought process as a way to show difficulties in designing the extension, and hopefully uncover where the configuration system/scheme/patterns can be improved. 

### Initial graph 
![image](mp_splat.png)

The above image shows the relevant portion of the graph I am planning to extend. Further to the right in the graph (i.e. after the tensor packer), the available tensors get written to a file. These tensors are the MP2 and MP3 tensors after going through CellViews of a resampled & thresholded splat frame. 

### Description of what I want to achieve
I would like to include the resampled & thresholded splat frame in the output. At some point, I need to concatenate the data along the channel dimension (like I do for the CellView'd data). I can't do this before putting it into cell views, so I need to fan the data out at some point. We have 2 clear options where this can occur:

### Graph 'Strategies'

#### 1: Fanout TensorSet immediately before CellViews (after packing)
![image](mp_splat_option1.png)

#### 2: Fanout Tensors before packing
![image](mp_splat_option2.png)

There isn't any clear winner yet. Option 2 is appealing because we don't have to unpack before concatenating, though I require 2 more fanouts. I don't think there's any real practical difference in terms of performance etc. I haven't looked too deeply into the configuration with this all in mind yet, but at this point, I think option 1 would be better because I don't have to worry about the fanouts. We'll see whether that's true. 

### Analyzing configuration
Below is the relevant portion of the config. I resample each group (U, V, W0, W1) then reduce W0&W1 -> Shared W view. So the shuntlines is rank 3 (views).
<pre>
    /// Put some of the previously-defined things together 
    /// This is one sticking point in the UI schema that I don't really like or at least
    /// took me a second to wrap my head around. What we're is defining progressively larger
    /// portions of the subgraph. The way we define that is by choosing some set of edges
    /// in the existing graph that we then have to splice and graft into in a way.
    threshold_cellviews_for_splat(out_views=[0,1,2], chunk_size=0, extra_name="_splat"):: 
        local this_name = $.this_name(extra_name, '_threshold_cellviews');

        local cellviews = $.cellviews_tensors(out_views=out_views, chunk_size=chunk_size, extra_name=extra_name);
        
        local name = $.this_name(extra_name, "_split");

        local stack = pg.pnode({
            type: 'SPNGReduce',
            name: $.this_name(extra_name, '_stack'),
            data: {
                multiplicity: 3,
                operation: "cat",
                dim: -2,
            } + control,
        }, nin=3, nout=1);

        pg.shuntlines([
            $.resample_group_to_views(extra_name=extra_name),
            pg.crossline($.thresholds_for_splat(extra_name=extra_name)),
            cellviews,
            stack,
        ]),
</pre>
Option 1 seeks to fanout the TensorSet immediately before the CellViews node. However, when I wrote this, I thought I'd save some space and take the 3 view tensors as input to the 'cellviews_tensors' function (it handles packing -> cellviews -> unpacking for me). Turns out, it would actually be nice to access the Packed TensorSet and fan it out. I could do that, and write packing/unpacking nodes explicitly here.
