# #for lax and guassian

serial_files = ['/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/serial/lax/2_threads_static,.txt']
for file in serial_files:
    with open(file, 'r') as curr_file:
        data = curr_file.readlines()
        guassian = data[8].split(": ")[1]
        print(guassian)
        lax_serial_time = data[10].split(": ")[1]
        print(lax_serial_time)
        with open("guassian_logs.txt", 'a') as file_to_write:
            file_to_write.write(guassian)
            file_to_write.write("\n")
        with open("lax_logs.txt", 'a') as file_to_write:
            file_to_write.write(lax_serial_time)
            file_to_write.write("\n")

files = ['/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/2_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/4_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/6_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/8_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/12_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/10_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/16_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/20_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/32_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/45_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/lax/64_threads_static.txt']
for file in files:
    with open(file, 'r') as current_file:
        data = current_file.readlines()
        guassian = data[8].split(": ")[1]
        print(guassian)
        lax_parallel_time = data[11].split(": ")[1]
        print(lax_parallel_time)
        with open("guassian_logs.txt", 'a') as file_to_write:
            file_to_write.write(guassian)
            file_to_write.write("\n")
        with open("lax_logs.txt", 'a') as file_to_write:
            file_to_write.write(lax_parallel_time)
            file_to_write.write("\n")
        
#for first_order
serial_files = ['/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/serial/first_order/2_threads_static,.txt']
for file in serial_files:
    with open(file, 'r') as curr_file:
        data = curr_file.readlines()
        serial_time = data[10].split(": ")[1]
        print(serial_time)
        with open("first_order_logs.txt", 'a') as file_to_write:
            file_to_write.write(serial_time)
            file_to_write.write("\n")

files = ['/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/2_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/4_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/6_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/8_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/10_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/12_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/16_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/20_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/32_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/45_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/first_order/64_threads_static.txt']
for file in files:
    with open(file, 'r') as current_file:
        print(file)
        data = current_file.readlines()
        parallel_time = data[11].split(": ")[1]
        print(parallel_time)
        with open("first_order_logs.txt", 'a') as file_to_write:
            file_to_write.write(parallel_time)
            file_to_write.write("\n")
            
#for second_order
serial_files = ['/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/serial/second_order/2_threads_static,.txt']
for file in serial_files:
    with open(file, 'r') as curr_file:
        data = curr_file.readlines()
        serial_time = data[10].split(": ")[1]
        with open("second_order_logs.txt", 'a') as file_to_write:
            file_to_write.write(serial_time)
            file_to_write.write("\n")

files = ['/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/2_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/4_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/6_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/8_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/10_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/12_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/16_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/20_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/32_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/45_threads_static,.txt',
         '/home/toberoi/hpc/project-1-winter-2023-tinaoberoi/logs/200N/parallel/second_order/64_threads_static.txt']
for file in files:
    with open(file, 'r') as current_file:
        data = current_file.readlines()
        parallel_time = data[11].split(": ")[1]
        print(parallel_time)
        with open("second_order_logs.txt", 'a') as file_to_write:
            file_to_write.write(parallel_time)
            file_to_write.write("\n")