writeline(file, s) = write(file, string(s, "\r\n"))
for i = 1:30
	filename = "ambulances_$i.csv"
	open(filename, "w") do file
		writeline(file, "ambulances,,")
		writeline(file, "index,stationIndex,class")
		for j = 1:i
			writeline(file, "$j,1,1")
		end
	end
end
