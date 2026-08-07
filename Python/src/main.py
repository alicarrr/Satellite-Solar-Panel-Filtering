import lcm 
from include.simMessages import sim_data_t

def sim_data_handler(channel, data):
    msg = sim_data_t.decode(data)
    print(f"Received data on channel: {channel}")
    print(msg.t)
    print(f"pos = {msg.pos}")


lc = lcm.LCM()
subscriber = lc.subscribe("TEST", sim_data_handler)

try:
    while True:
        lc.handle()
except KeyboardInterrupt:
    pass