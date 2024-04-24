import ray
@ray.remote
def f(x):
    print('coucou2')
    return x * x

if __name__ == "__main__":


    ray.init()
    print('coucouc1')
    futures = [f.remote(i) for i in range(4)]
    print(ray.get(futures)) # [0, 1, 4, 9]
    ray.close()
