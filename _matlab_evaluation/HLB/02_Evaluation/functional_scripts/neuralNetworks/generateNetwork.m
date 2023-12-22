function net = generateNetwork(numChannels)
    % generating using dlnetwork
    layers = [
    sequenceInputLayer(numChannels)
    lstmLayer(128)
    fullyConnectedLayer(numChannels)
    regressionLayer];

    net = dlnetwork(layers);
end

