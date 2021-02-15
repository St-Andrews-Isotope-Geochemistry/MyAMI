classdef MyAMI < handle
    properties
        precalculated = NaN
        function_handles
        results
        use_cache
        cache = containers.Map();
        k_order = ["k0","k1","k2","kb","kw","kc","ka","ks"];
        method
    end
    methods
        function self = MyAMI(method,use_cache,previous_cache)
            if strcmp(method,"MyAMI")
                self.method = "MyAMI";
            elseif strcmp(method,"Precalculated");
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
                    self.precalculated = readmatrix(self.getMyAMIPath()+"Precalculated.xls","Range","A:BV");
                end
                if strcmp(self.method,"MyAMI")
                    [self.results,~] = self.run(temperature,salinity,calcium,magnesium);
                elseif strcmp(self.method,"Precalculated")
                    self.results = self.fromPrecalculated(temperature,salinity,calcium,magnesium,self.precalculated,self.function_handles);
                end
                if self.use_cache
                    self.cache(self.getKey(temperature,salinity,calcium,magnesium)) = self.results;
                end
            end
        end
        function clear(self)
            self.results = NaN;
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
%             key = strrep(initial_key,".","");
        end
        function [k_values,k_values_correction_output] = run(temperature,salinity,calcium,magnesium)
            command = join(["python",MyAMI.MyAMI.getMyAMIPath()+"/PITZER.py ",temperature,salinity,calcium,magnesium]," ");
            [status,result] = system(command);
            if status==0
            % Should check status and result for reasonable values
                result_cleaned = erase(string(result(2:end-2)),newline);
                result_values = str2num(result_cleaned);
                input_values = result_values(1:4);
                output_values = result_values(5:end-1); % Last value unused in the examples...
                output_matrix = reshape(output_values,[],3);
                
                k_order = ["k0","k1","k2","kb","kw","kc","ka","ks"];
                k_values_correction = output_matrix(:,1)./output_matrix(:,3);
                
                k_values_correction_output = containers.Map(k_order,output_matrix(:,1)./output_matrix(:,3));
                k_values = containers.Map(k_order,output_matrix(:,2).*k_values_correction);
            else
                error(join(["Problem running Python:",result]," "));
            end
        end
        function [k_values] = fromPrecalculated(temperature,salinity,calcium,magnesium,precalculated,functions)
            if isstring(precalculated)
                precalculated = readmatrix(precalculated,"Range","A:BV");
            end
            
            ionic_strength = (19.924*salinity)/(1000-1.005*salinity);
            
            calcium_magnesium = precalculated(:,1:2);
            header_rows = sum(isnan(calcium_magnesium(:,1)));
            
            calcium_known = calcium_magnesium(header_rows+1:end,1);
            magnesium_known = calcium_magnesium(header_rows+1:end,2);
            
            calcium_unique = unique(calcium_known);
            magnesium_unique = unique(magnesium_known);
            
            calcium_resolution = calcium_unique(2)-calcium_unique(1);
            calcium_range = max(calcium_unique)-min(calcium_unique);
            calcium_minimum = min(calcium_unique);
            
            magnesium_resolution = magnesium_unique(2)-magnesium_unique(1);
            magnesium_range = max(magnesium_unique)-min(magnesium_unique);
            magnesium_minimum = min(magnesium_unique);
            
            if mod(calcium,calcium_resolution)==0 && mod(magnesium,magnesium_resolution)==0 % There's an exact result in the spreadsheet
                index = header_rows+(1+magnesium/magnesium_resolution)+((calcium/calcium_resolution)*(1+magnesium_range/magnesium_resolution));
                coefficients = [precalculated(round(index),4:74),NaN];
            else
                if mod(calcium,calcium_resolution)==0
                    calcium_query = [floor(calcium*1000)/1000,ceil((calcium+1e-6)*1000)/1000];
                    magnesium_query = [floor(magnesium*1000)/1000,ceil(magnesium*1000)/1000];
                elseif mod(magnesium,magnesium_resolution)==0
                    calcium_query = [floor(calcium*1000)/1000,ceil(calcium*1000)/1000];
                    magnesium_query = [floor(magnesium*1000)/1000,ceil((magnesium+1e-6)*1000)/1000];
                else
                    calcium_query = [floor(calcium*1000)/1000,ceil(calcium*1000)/1000];
                    magnesium_query = [floor(magnesium*1000)/1000,ceil(magnesium*1000)/1000];
                end
                
                for calcium_query_index = 1:numel(calcium_query)
                    for magnesium_query_index = 1:numel(magnesium_query)
                        index = round(header_rows+(1+(magnesium_query(magnesium_query_index)-magnesium_minimum)/magnesium_resolution)+(((calcium_query(calcium_query_index)-calcium_minimum)/calcium_resolution)*(1+magnesium_range/magnesium_resolution)));
                        raw_coefficients(magnesium_query_index,calcium_query_index,:) = precalculated(index,4:74);
                    end
                end
                
                for coefficient_index = 1:size(raw_coefficients,3)
                    interpolated_coefficients(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * raw_coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                end
                coefficients = [interpolated_coefficients,NaN];
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
end