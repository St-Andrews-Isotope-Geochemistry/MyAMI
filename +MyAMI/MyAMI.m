classdef MyAMI < handle
    properties
        precalculated = NaN
        function_handles
        results
        use_cache
        cache = containers.Map();
        k_order = ["k0","k1","k2","kb","kw","kc","ka","ks"];
        method
        header_rows
        calcium_resolution
        calcium_minimum
        calcium_maximum
        magnesium_resolution
        magnesium_minimum
        magnesium_maximum
    end
    methods
        function self = MyAMI(method,use_cache,previous_cache)
            if strcmp(method,"MyAMI")
                self.method = "MyAMI";
            elseif strcmp(method,"Precalculated")
                self.method = "Precalculated";
            else
                error("Method unknown");
            end
            
            if nargin>1
                if use_cache==true
                    self.use_cache = true;
                elseif use_cache==false
                    self.use_cache = false;
                else
                    error("Cache must be true or false");
                end
            else
                self.use_cache = true;
            end
            
            if nargin>2
                self.cache = previous_cache;
            end
            
            k0_function = @(coefficients,t,s,i) coefficients(1) + (100*coefficients(2))/t + coefficients(3)*log(t/100) + s*(coefficients(4) + (coefficients(5)*t)/100 + coefficients(6)*(t/100)^2);
            k1_function = @(coefficients,t,s,i) log(10^(coefficients(1) + coefficients(2)/t + coefficients(3)*log(t) + coefficients(4)*s + coefficients(5)*s^2));
            k2_function = @(coefficients,t,s,i) log(10^(coefficients(1) + coefficients(2)/t + coefficients(3)*log(t) + coefficients(4)*s + coefficients(5)*s^2));
            kb_function = @(coefficients,t,s,i) coefficients(1) + coefficients(2)*sqrt(s) + coefficients(3)*s + (1/t)*(coefficients(4) + coefficients(5)*sqrt(s) + coefficients(6)*s + coefficients(7)*s^1.5 + coefficients(8)*s^2) + log(t)*(coefficients(9) + coefficients(10)*sqrt(s) + coefficients(11)*s) + coefficients(12)*t*sqrt(s);
            kw_function = @(coefficients,t,s,i) coefficients(1) + coefficients(2)/t + coefficients(3)*log(t) + sqrt(s)*(coefficients(4)/t + coefficients(5) + coefficients(6)*log(t)) + coefficients(7)*s;
            kc_function = @(coefficients,t,s,i) log(10^(coefficients(1) + coefficients(2)*t + coefficients(3)/t + coefficients(4)*log10(t) + sqrt(s)*(coefficients(5) + coefficients(6)*t + coefficients(7)/t) + coefficients(8)*s + coefficients(9)*s^1.5));
            ka_function = @(coefficients,t,s,i) log(10^(coefficients(1) + coefficients(2)*t + coefficients(3)/t + coefficients(4)*log10(t) + sqrt(s)*(coefficients(5) + coefficients(6)*t + coefficients(7)/t) + coefficients(8)*s + coefficients(9)*s^1.5));
            ks_function = @(coefficients,t,s,i) coefficients(1) + coefficients(2)/t + coefficients(3)*log(t) + sqrt(i)*(coefficients(4)/t + coefficients(5) + coefficients(6)*log(t)) + i*(coefficients(7)/t + coefficients(8) + coefficients(9)*log(t)) + (coefficients(10)*i^1.5)/t + (coefficients(11)*i^2)/t + log(1-0.001005*s);
            
            k_functions = {k0_function,k1_function,k2_function,kb_function,kw_function,kc_function,ka_function,ks_function};
            
            self.function_handles = containers.Map(self.k_order,k_functions);
        end
        function calculate(self,temperature,salinity,calcium,magnesium)
            if self.use_cache && isKey(self.cache,self.getKey(temperature,salinity,calcium,magnesium))
                self.results = self.cache(self.getKey(temperature,salinity,calcium,magnesium));
            else
                if isnan(self.precalculated)
                    self.precalculated = self.parsePrecalculated(self.getPrecalculated());
                end
                if strcmp(self.method,"MyAMI")
                    [self.results,~] = self.run(temperature,salinity,calcium,magnesium);
                elseif strcmp(self.method,"Precalculated")
                    self.results = self.fromPrecalculated(temperature,salinity,calcium,magnesium,self.precalculated,self.function_handles,[self.calcium_minimum,self.calcium_maximum,self.calcium_resolution],[self.magnesium_minimum,self.magnesium_maximum,self.magnesium_resolution]);
                end
                if self.use_cache
                    self.cache(self.getKey(temperature,salinity,calcium,magnesium)) = self.results;
                end
            end
        end
        function clear(self)
            self.results = NaN;
        end
        function precalculated = getPrecalculated(self)
            precalculated = readmatrix(self.getMyAMIPath()+"Precalculated.xls","Range","A:BV");
        end
        function precalculated_reshaped = parsePrecalculated(self,precalculated)
            if nargin<2 || (numel(precalculated)==1 && isnan(precalculated))
                precalculated = self.getPrecalculated();
            end
            calcium_magnesium = precalculated(:,1:2);
            self.header_rows = sum(isnan(calcium_magnesium(:,1)));
            
            calcium_known = calcium_magnesium(self.header_rows+1:end,1);
            magnesium_known = calcium_magnesium(self.header_rows+1:end,2);
            
            calcium_unique = unique(calcium_known);
            magnesium_unique = unique(magnesium_known);
            
            self.calcium_resolution = calcium_unique(2)-calcium_unique(1);
            self.calcium_minimum = min(calcium_unique);
            self.calcium_maximum = max(calcium_unique);
            
            self.magnesium_resolution = magnesium_unique(2)-magnesium_unique(1);
            self.magnesium_minimum = min(magnesium_unique);
            self.magnesium_maximum = max(magnesium_unique);
            
            precalculated_reshaped = permute(reshape(precalculated(self.header_rows+1:end,:),numel(calcium_unique),numel(magnesium_unique),[]),[2,1,3]);
        end
        function [k_values,k_values_correction_output] = run(self,temperature,salinity,calcium,magnesium)
            command = join(["python",MyAMI.MyAMI.getMyAMIPath()+"/PITZER.py ",temperature,salinity,calcium,magnesium]," ");
            [status,result] = system(command);
            if status==0
            % Should check status and result for reasonable values
                result_cleaned = erase(string(result(2:end-2)),newline);
                result_values = str2num(result_cleaned);
                input_values = result_values(1:4);
                output_values = result_values(5:end-1); % Last value unused in the examples...
                output_matrix = reshape(output_values,[],3);
                
                k_order = ["kw","k1","k2","kc","kb","ka","k0","ks"];
                k_values_correction = output_matrix(:,1)./output_matrix(:,3);
                
                k_values_correction_output = containers.Map(k_order,output_matrix(:,1)./output_matrix(:,3));
                k_values = containers.Map(k_order,output_matrix(:,2).*k_values_correction);
            else
                error(join(["Problem running Python:",result]," "));
            end
        end
        function [k_values] = fromPrecalculated(self,temperature,salinity,calcium,magnesium,precalculated,functions,calcium_information,magnesium_information)
            if isstring(precalculated)
                raw_precalculated = readmatrix(precalculated,"Range","A:BV");
            
                calcium_magnesium = raw_precalculated(:,1:2);
                self.header_rows = sum(isnan(calcium_magnesium(:,1)));
                
                calcium_known = calcium_magnesium(self.header_rows+1:end,1);
                magnesium_known = calcium_magnesium(self.header_rows+1:end,2);
                
                calcium_unique = unique(calcium_known);
                magnesium_unique = unique(magnesium_known);
                
                calcium_minimum = min(calcium_unique);
                calcium_maximum = max(calcium_unique);
                calcium_resolution = calcium_unique(2)-calcium_unique(1);
                
                magnesium_minimum = min(magnesium_unique);
                magnesium_maximum = max(magnesium_unique);
                magnesium_resolution = magnesium_unique(2)-magnesium_unique(1);
                
                precalculated = reshape(raw_precalculated(self.header_rows+1:end,:),numel(calcium_unique),numel(magnesium_unique),[]);
            else
                calcium_minimum = calcium_information(1);
                calcium_maximum = calcium_information(2);
                calcium_resolution = calcium_information(3);
                
                magnesium_minimum = magnesium_information(1);
                magnesium_maximum =  magnesium_information(2);
                magnesium_resolution = magnesium_information(3);
            end
            
            ionic_strength = (19.924*salinity)/(1000-1.005*salinity);
            
            if calcium>=calcium_minimum && calcium<=calcium_maximum && magnesium>=magnesium_minimum && magnesium<=magnesium_maximum % Interpolate
                if mod(calcium,calcium_resolution)==0 && mod(magnesium,magnesium_resolution)==0 % There's an exact result in the spreadsheet
                    index = round([1+((calcium-calcium_minimum)/calcium_resolution),1+((magnesium-magnesium_minimum)/magnesium_resolution)]);
                    coefficients = [squeeze(precalculated(index(1),index(2),4:74))',NaN];
                else % No exact match
                    if mod(calcium,calcium_resolution)==0 % Calcium match but no magnesium match
                        magnesium_query = [floor(magnesium*1000)/1000,ceil(magnesium*1000)/1000];
                        
                        raw_coefficients = NaN(1,2,1+74-4);
                        for magnesium_query_index = 1:numel(magnesium_query)
                            index = round([1+((calcium-calcium_minimum)/calcium_resolution),1+((magnesium_query(magnesium_query_index)-magnesium_minimum)/magnesium_resolution)]);
                            raw_coefficients(1,magnesium_query_index,:) = precalculated(index(1),index(2),4:74);
                        end                        
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = self.linearEstimate(magnesium_query,raw_coefficients(:,:,coefficient_index),magnesium);
                        end
                    elseif mod(magnesium,magnesium_resolution)==0 % Magnesium match but not calcium match
                        calcium_query = [floor(calcium*1000)/1000,ceil(calcium*1000)/1000];
                        
                        raw_coefficients = NaN(2,1,1+74-4);
                        for calcium_query_index = 1:numel(calcium_query)
                            index = round([1+((calcium_query(calcium_query_index)-calcium_minimum)/calcium_resolution),1+((magnesium-magnesium_minimum)/magnesium_resolution)]);
                            raw_coefficients(calcium_query_index,1,:) = precalculated(index(1),index(2),4:74);
                        end
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = self.linearEstimate(calcium_query,raw_coefficients(:,:,coefficient_index)',calcium);
                        end                        
                    else
                        calcium_query = [floor(calcium*1000)/1000,ceil(calcium*1000)/1000];
                        magnesium_query = [floor(magnesium*1000)/1000,ceil(magnesium*1000)/1000];
                        
                        raw_coefficients = NaN(2,2,1+74-4);
                        for calcium_query_index = 1:numel(calcium_query)
                            for magnesium_query_index = 1:numel(magnesium_query)
                                index = round([1+((calcium_query(calcium_query_index)-calcium_minimum)/calcium_resolution),1+((magnesium_query(magnesium_query_index)-magnesium_minimum)/magnesium_resolution)]);
                                raw_coefficients(calcium_query_index,magnesium_query_index,:) = precalculated(index(1),index(2),4:74);
                            end
                        end
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * raw_coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                        end
                    end
                    coefficients = interpolated_coefficients;
                end
            else % Extrapolate
                if calcium>=calcium_minimum && calcium<=calcium_maximum
                    magnesium_query = [magnesium_maximum-magnesium_resolution,magnesium_maximum];
                    if mod(calcium,calcium_resolution)==0
                        for calcium_query_index = 1:numel(calcium_query)
                            for magnesium_query_index = 1:numel(magnesium_query)
                                index = round([1+((calcium-calcium_minimum)/calcium_resolution),1+((magnesium_query(magnesium_query_index)-magnesium_minimum)/magnesium_resolution)]);
                                raw_coefficients(calcium_query_index,magnesium_query_index,:) = precalculated(index(1),index(2),4:74);
                            end
                        end
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * raw_coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                        end
                    else
                        calcium_query = [floor(calcium*1000)/1000,ceil(calcium*1000)/1000];
                        for calcium_query_index = 1:numel(calcium_query)
                            for magnesium_query_index = 1:numel(magnesium_query)
                                index = round([1+((calcium_query(calcium_query_index)-calcium_minimum)/calcium_resolution),1+((magnesium_query(magnesium_query_index)-magnesium_minimum)/magnesium_resolution)]);
                                raw_coefficients(calcium_query_index,magnesium_query_index,:) = precalculated(index(1),index(2),4:74);
                            end
                        end
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * raw_coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                        end
                    end
                elseif magnesium>=magnesium_minimum && magnesium<=magnesium_maximum
                    calcium_query = [calcium_maximum-calcium_resolution,calcium_maximum];
                    if mod(magnesium,magnesium_resolution)==0
                        for calcium_query_index = 1:numel(calcium_query)
                            index = round([1+((calcium_query(calcium_query_index)-calcium_minimum)/calcium_resolution),1+((magnesium-magnesium_minimum)/magnesium_resolution)]);
                            raw_coefficients(calcium_query_index,magnesium_query_index,:) = precalculated(index(1),index(2),4:74);
                        end
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * raw_coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                        end
                    else
                        magnesium_query = [floor(magnesium*1000)/1000,ceil(magnesium*1000)/1000];
                        for calcium_query_index = 1:numel(calcium_query)
                            for magnesium_query_index = 1:numel(magnesium_query)
                                index = round([1+((calcium_query(calcium_query_index)-calcium_minimum)/calcium_resolution),1+((magnesium_query(magnesium_query_index)-magnesium_minimum)/magnesium_resolution)]);
                                raw_coefficients(calcium_query_index,magnesium_query_index,:) = precalculated(index(1),index(2),4:74);
                            end
                        end
                        
                        interpolated_coefficients = NaN(2+74-4,1);
                        for coefficient_index = 1:size(raw_coefficients,3)
                            interpolated_coefficients(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * raw_coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                        end
                    end
                end
                coefficients = interpolated_coefficients;
            end
            coefficient_map = containers.Map();
            k_order = ["k0","k1","k2","kb","kw","kc","ka","ks"];
            k_count = 1;
            coefficient_start = 1;
            for coefficient_index = 1:numel(coefficients)
                if ~isnan(coefficients(coefficient_index))
                    continue
                else
                    coefficient_map(k_order(k_count)) = coefficients(coefficient_start:coefficient_index-1);
                    k_count = k_count+1;
                    coefficient_start = coefficient_index+1;
                end
            end
            
            output_map = containers.Map();
            for ck = k_order
                current_function = functions(ck);
                output_map(ck) = exp(current_function(coefficient_map(ck),temperature+273.15,salinity,ionic_strength));
            end
            k_values = output_map;
        end
        
    end
    methods (Static=true)
        function MyAMI_path = getMyAMIPath()
            MyAMI_search = what("+MyAMI");
            current_directory = pwd; 
            MyAMI_path = strrep(strrep(strrep(join([MyAMI_search(1).path],""),current_directory,"."),"\","/"),"+MyAMI","");
        end
        function key = getKey(temperature,salinity,calcium,magnesium)
            key = string(temperature)+string(salinity)+string(calcium)+string(magnesium);
        end
        function estimate = linearEstimate(boundaries,to_interpolate,value)
            assert(numel(boundaries)==2,"Must be 2 boundaries");
            assert(size(to_interpolate,2)==2,"Matrix to interpolate must be n x 2 x m");
            assert(numel(size(to_interpolate))<=2,"Matrix to interpolate must be 1D or 2D");
            estimate = (1/(boundaries(2)-boundaries(1))).*([boundaries(2)-value,value-boundaries(1)]*to_interpolate')';
        end
    end
end