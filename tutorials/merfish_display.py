import intake
import param
import panel as pn
import holoviews as hv
import panel.widgets as pnw
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr
import zarr
from pathlib import Path
from holoviews.operation.datashader import regrid
from imaxt_image.image import zscale



class Catalogue(param.Parameterized):
    
    gene = param.ObjectSelector()
    exact = param.Boolean(default=True, doc='Exact bitcode')
    
    def __init__(self, catalogue):
        self.ds = intake.open_catalog(catalogue)
        self.df = self.ds.merged.read()
        
        self._genes = sorted(self.df.gene.unique().tolist())
        self.param['gene'].objects=self._genes
        self.param['gene'].default=self._genes[0]
        super().__init__()
        
    def head(self):
        return self.df.head()
    
    def _get_spatial_distribution(self, gene, exact):
        data = self.df[self.df['gene']==gene]
        if exact:
            data = data[data.distance==0]
        points = hv.Scatter((data['x_abs'], -data['y_abs']), 'x_abs', 'y_abs')
        res = points.opts(size=2, title=self.gene)
        return res # main = res.opts(width=800, height=400, title=self.gene)
        
    def _get_fpkm(self, gene, exact):
        data = self.df
        if exact:
            data = data[data.distance==0]
        csv = pd.read_csv('4t1_tissue_rnaseq_abundance.csv')
        e = data.groupby('gene').size().reset_index(name='n')
        res = csv.join(e.set_index('gene'), on='gene_name')
        fpkm = hv.Scatter(data=res, kdims=['FPKM',], vdims=['n',]).opts(logx=True, logy=True) * \
               hv.Scatter(data=res[res.gene_name.str.contains(gene)], kdims=['FPKM',], vdims=['n',]).opts(size=10, color='blue')

        r = res[res.FPKM>0]
        pnumber = pearsonr(np.log10(r.FPKM), np.log10(r.n))
        fpkm = fpkm.opts(title=f'r: {pnumber[0]:.2f} p: {pnumber[1]:.2g}')
        
        #hover = HoverTool(tooltips=[
        #('desc', '@desc'),
        #])
        #source = ColumnDataSource(data=dict(x=xx, y=yy, desc=desc))

        #p = figure(title=f'r: {pnumber[0]:.2f} p: {pnumber[1]:.2g}', plot_width=400, plot_height=400, y_axis_type='log', x_axis_type='log', tools=[hover])
        #p.circle('x', 'y', source=source)
        return fpkm

            
    @param.depends('gene', 'exact')
    def view(self):
        spatial = self._get_spatial_distribution(self.gene, self.exact)
        fpkm = self._get_fpkm(self.gene, self.exact)
        expression = self._get_expression(self.gene, self.exact)
        
        gspec = pn.GridSpec(width=1000, height=800)
        gspec[0, 0:6] = spatial
        gspec[1, 0:2] = fpkm
        gspec[1, 2:6] = expression
        
        return gspec
        
    def _get_expression(self, gene, exact):    
        data = self.df
        if self.exact:
            data = data[data.distance==0]
        data = data.groupby('gene').size().reset_index(name='n')
        data = data.sort_values('n', ascending=False)
        bb = hv.Bars(data, kdims=('gene'), vdims=('n')).opts(line_color='steelblue', fill_color='steelblue', bar_width=1, xaxis=None, logy=True) * \
        hv.Bars(data[data.gene.str.contains('blank')], kdims=('gene'), vdims=('n')).opts(line_color='red', fill_color='red', bar_width=1, xaxis=None, logy=True) * \
        hv.Bars(data[data.gene.str.contains(gene)], kdims=('gene'), vdims=('n')).opts(line_color='red', bar_width=1, xaxis=None, logy=True)        
        return bb
    
    def display(self, height=400):
        res = pn.Column(self.param, self.view)
        return res
        
class ImageViewer:
    def __init__(self, path):
        self.path = Path(path)
        self.arr = zarr.open(path + '/zarr')
        
    def load_image(self, plane, channel, cycle, fov):
        g = [*self.arr[f'fov={fov}/z={plane}/cycle={cycle}'].groups()]
        g.sort()
        im1 = g[channel][1]['raw'][:]
        vmin, vmax = zscale(im1)
        im1 = im1.clip(vmin, vmax)
        im1 = hv.Image(im1, bounds=(0, 0, 2048, 2048), vdims='Intensity')
        return im1

    def display(self):
        # Define DynamicMap with z-dimension to slide through
        image_stack = hv.DynamicMap(self.load_image, kdims=['plane', 'channel', 'cycle', 'fov'])


        # Apply regridding in case data is large and set a global Intensity range 
        regridded = regrid(image_stack).redim.range(plane=(0,self.arr.attrs['planes']-1), channel=(0,2), cycle=(0,self.arr.attrs['cycles']-1), fov=(0,self.arr.attrs['fov']-1)) #+ hist_stack.redim.range(z=(0,3), c=(0,3))
        display_obj = regridded.opts(plot={'Image': dict(width=700, height=700, colorbar=False, 
                                                         xaxis=None, yaxis=None, aspect='equal')})
        return display_obj