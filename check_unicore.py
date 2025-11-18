try:
    import unicore
    print('Uni-Core 已安装')
except ImportError as e:
    print(f'Uni-Core 未安装: {e}')
