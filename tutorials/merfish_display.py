import intake
import param
import panel as pn
import holoviews as hv
import panel.widgets as pnw
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr


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
        