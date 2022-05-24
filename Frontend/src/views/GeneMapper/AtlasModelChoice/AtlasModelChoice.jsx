import HelpIcon from '@mui/icons-material/Help';
import { Grid, Typography, Stack, Alert, Box, Container, Tooltip } from "@mui/material";
import AtlasCardSelect from "components/Cards/AtlasCardSelect";
import { ModelCardSelect } from "components/Cards/ModelCardSelect";
import CustomButton from 'components/CustomButton';
import { useHistory} from 'react-router-dom';
import { useEffect, useState } from 'react';
import ProjectMock from "shared/services/mock/projects"

import Clear from '@mui/icons-material/Clear';
import ArrowForwardIcon from '@mui/icons-material/ArrowForward';

import ModelService from 'shared/services/Model.service';
import AtlasService from 'shared/services/Atlas.service';

import { colors } from 'shared/theme/colors';

function AtlasModelChoice({ 
    activeStep, setActiveStep, 
    selectedAtlas, setSelectedAtlas, 
    selectedModel, setSelectedModel, steps, path,
    compatibleModels
}) {
    let [atlases, setAtlases] = useState([]);
    let [models, setModels] = useState([]);
    const [showWarning, setShowWarning] = useState(false);
    const history = useHistory();

    let headerStyle = {
        color: "#003560",
        fontSize: "1.6rem",
        fontWeight: "bold",
    }
    
    useEffect(() => {
        AtlasService.getAtlases().then(a => {
            a.map(a => {
                //adjust numberOfCells
                let numberOfCells = a.numberOfCells;
                let dimension = ""
                if (numberOfCells > 1000000000) {
                    numberOfCells = Math.round(numberOfCells / 1000000000);
                    dimension = "B";
                } else if (numberOfCells > 1000000) {
                    numberOfCells = Math.round(numberOfCells / 1000000);
                    dimension = "M";
                } else if (numberOfCells > 1000) {
                    numberOfCells = Math.round(numberOfCells / 1000);
                    dimension = "K"
                }
                a.numberOfCells = numberOfCells + dimension;

                //adjust modalities
                if (!(typeof a.modalities === "string")) {
                    //modalities ist array of strings
                    if (a.modalities.length == 0) {
                        a.modalities = "None";
                    } else if (a.modalities.length == 1) {
                        a.modalities = a.modalities[0];
                    } else {
                        a.modalities = a.modalities[0] + ', ...';
                    }
                }
                return a;
            })
            setAtlases(a);
        });
        ModelService.getModels().then(m => {
            m.map(model => {
                model.requirements = ['Batch/study information is mandatory and should be labeled as “batch”'];
                if (model.name === 'scVI') {
                    model.requirements.push('Cell type information should be labeled as “cell_type”');
                    model.requirements.push('For unlabeled cells, the value for “cell_type” should be “Unknown”');
                }
                if (model.name === 'scANVI') {
                    model.requirements.push('Cell type information should be labeled as “cell_type”');
                }
            })
            setModels(m);
        });
    }, []);

    return (
        <div>
            {showWarning &&
            <Alert severity="error">
                Select an Atlas and a fitting Model before continuing
            </Alert>}
            <Typography 
            variant="h5" 
            sx={{ 
                fontWeight: 'bold', 
                pb: "1em", 
                mt: "1.5em",
            }}>
                Pick an Atlas 
            </Typography>

            <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
                {
                    atlases.map(a => 
                        <Grid item height="320px">
                            <AtlasCardSelect width="225px"
                                height="97%"
                                title={a.name}
                                modalities={a.modalities}
                                cellsInReference={a.numberOfCells}
                                species={a.species}
                                imgLink={a.previewPictureURL}
                                selected={selectedAtlas.name===a.name}
                                onSelect={setSelectedAtlas}
                                atlasObject={a}
                            />      
                        </Grid>
                    )
                }
            </Grid>
            
            <Typography variant="h5" sx={{ fontWeight: 'bold', pb: "1em" }} marginTop="32px">
                Pick a Model
            </Typography>

            <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
            {
                models.map(m => 
                    <Grid item height="150px">
                    <ModelCardSelect width="225px"
                            height="97%"
                            title={m.name} 
                            description={m.description}
                            selected={selectedModel.name===m.name}
                            onSelect={setSelectedModel}
                            modelObject={m}
                            disabled={!compatibleModels.map(m => m.toLowerCase()).includes(m.name.toLowerCase()) || compatibleModels.length == 0}
                    />               
                </Grid>)
            }
            </Grid>
            <Stack direction='row' justifyContent='space-between' sx={{ marginTop:'20px', marginBottom: "3em"}}>
                <CustomButton type='tertiary' onClick={() => history.push(`${path}`)}>
                <Clear/>&nbsp; Cancel
                </CustomButton>
                <Tooltip title={!selectedAtlas || !selectedModel ? "Select an Atlas and a fitting Model before continuing" : ''} placement='top'>
                    <Box
                        onClick={!selectedAtlas || !selectedModel ? () => setShowWarning(true) : ()=>{}}
                    >
                        <CustomButton 
                            type='primary' 
                            disabled={!selectedAtlas || !selectedModel} 
                            onClick={() => setActiveStep(1)}
                        >
                        Confirm&nbsp;<ArrowForwardIcon/>
                        </CustomButton>
                    </Box>
                </Tooltip>
            </Stack>
        </div>
    )
}

export default AtlasModelChoice