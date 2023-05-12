import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.colors import ListedColormap

import base64
from io import BytesIO


def get_graph():
    buffer = BytesIO()
    plt.savefig(buffer, format = 'svg')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()
    return graph

def volcano(df="dataframe", lfc=None, pv=None, lfc_thr=(1, 1), pv_thr=(0.05, 0.05), color=("green", "grey", "red"),
                valpha=1, geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False, figtype='png', axtickfontsize=9,
                axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None, plotlegend=False, legendpos='best',
                figname='volcano', legendanchor=None,
                legendlabels=['significant up', 'not significant', 'significant down'], theme=None):

    plt.switch_backend('AGG')

    _x = r'$ log_{2}(Fold Change)$'
    _y = r'$ -log_{10}(P-value)$'

    # plt.figure(figsize=(10, 20))

    color = color
    # check if dataframe contains any non-numeric character
    # assert general.check_for_nonnumeric(df[lfc]) == 0, 'dataframe contains non-numeric values in lfc column'
    # assert general.check_for_nonnumeric(df[pv]) == 0, 'dataframe contains non-numeric values in pv column'
    # this is important to check if color or logpv exists and drop them as if you run multiple times same command
    # it may update old instance of df
    df = df.drop(['color_add_axy', 'logpv_add_axy'], axis=1, errors='ignore')
    assert len(set(color)) == 3, 'unique color must be size of 3'
    df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr[0]), 'color_add_axy'] = color[0]  # upregulated
    df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr[1]), 'color_add_axy'] = color[2]  # downregulated
    df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate
    df['logpv_add_axy'] = -(np.log10(df[pv]))
    # plot
    assign_values = {col: i for i, col in enumerate(color)}
    color_result_num = [assign_values[i] for i in df['color_add_axy']]
    # assert len(set(color_result_num)) == 3, \
    #     'either significant or non-significant genes are missing; try to change lfc_thr or pv_thr to include ' \
    #     'both significant and non-significant genes'
    if theme == 'dark':
        general.dark_bg()
    plt.subplots(figsize=dim)
    if plotlegend:
        s = plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                        s=dotsize, marker=markerdot)
        assert len(legendlabels) == 3, 'legendlabels must be size of 3'
        plt.legend(handles=s.legend_elements()[0], labels=legendlabels, loc=legendpos, bbox_to_anchor=legendanchor)
    else:
        plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                    s=dotsize, marker=markerdot)
    if sign_line:
        plt.axhline(y=-np.log10(pv_thr[0]), linestyle='--', color='#7d7d7d', linewidth=1)
        plt.axvline(x=lfc_thr[0], linestyle='--', color='#7d7d7d', linewidth=1)
        plt.axvline(x=-lfc_thr[1], linestyle='--', color='#7d7d7d', linewidth=1)
    # GeneExpression.gene_plot(df, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle)

    if genenames is not None and genenames == "deg":
            for i in df[geneid].unique():
                if (df.loc[df[geneid] == i, lfc].iloc[0] >= lfc_thr[0] and df.loc[df[geneid] == i, pv].iloc[0] < pv_thr[0]) or \
                        (df.loc[df[geneid] == i, lfc].iloc[0] <= -lfc_thr[1] and df.loc[df[geneid] == i, pv].iloc[0] < pv_thr[1]):
                    if gstyle == 1:
                        plt.text(df.loc[df[geneid] == i, lfc].iloc[0], df.loc[df[geneid] == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(i, xy=(df.loc[df[geneid] == i, lfc].iloc[0], df.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
    elif genenames is not None and type(genenames) is tuple:
        for i in df[geneid].unique():
            if i in genenames:
                if gstyle == 1:
                    plt.text(df.loc[df[geneid] == i, lfc].iloc[0], df.loc[df[geneid] == i, 'logpv_add_axy'].iloc[0], i,
                                  fontsize=gfont)
                elif gstyle == 2:
                    plt.annotate(i, xy=(df.loc[df[geneid] == i, lfc].iloc[0], df.loc[df[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                 xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                 bbox=dict(boxstyle="round", alpha=0.1),
                                 arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                else:
                    print("Error: invalid gstyle choice")
                    sys.exit(1)
    elif genenames is not None and type(genenames) is dict:
        for i in df[geneid].unique():
            if i in genenames:
                if gstyle == 1:
                    plt.text(df.loc[df[geneid] == i, lfc].iloc[0], df.loc[df[geneid] == i, 'logpv_add_axy'].iloc[0],
                                  genenames[i], fontsize=gfont)
                elif gstyle == 2:
                    plt.annotate(genenames[i], xy=(df.loc[df[geneid] == i, lfc].iloc[0], df.loc[df[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                 xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                 bbox=dict(boxstyle="round", alpha=0.1),
                                 arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                else:
                    print("Error: invalid gstyle choice")
                    sys.exit(1)

    if axxlabel:
        _x = axxlabel
    if axylabel:
        _y = axylabel
    general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
    general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)


    volcano_plot = get_graph()
    return volcano_plot

    # class GeneExpression:

#     def __init__(self):
#         pass

#     @staticmethod
#     def gene_plot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle):
#
class general:
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')




    @staticmethod
    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        # plt.xticks(fontsize=9, fontname="sans-serif")
        # plt.yticks(fontsize=9, fontname="sans-serif")

    @staticmethod
    def axis_ticks(xlm=None, ylm=None, axtickfontsize=None, axtickfontname=None, ar=None):
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    @staticmethod
    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

    @staticmethod
    def check_for_nonnumeric(pd_series=None):
        if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
            return 0
        else:
            return 1

    @staticmethod
    def pvalue_symbol(pv=None, symbol=None):
        if 0.05 >= pv > 0.01:
            return symbol
        elif 0.01 >= pv > 0.001:
            return 2 * symbol
        elif pv <= 0.001:
            return 3 * symbol
        else:
            return None

    @staticmethod
    def get_file_from_gd(url=None):
        get_path = 'https://drive.google.com/uc?export=download&id=' + url.split('/')[-2]
        return pd.read_csv(get_path, comment='#')

    @staticmethod
    def dark_bg():
        plt.style.use('dark_background')


