import HelpIcon from '@mui/icons-material/Help';
import { Grid, Typography, Stack } from "@mui/material";
import AtlasCardSelect from "components/Cards/AtlasCardSelect";
import { ModelCardSelect } from "components/Cards/ModelCardSelect";
import CustomButton from 'components/CustomButton';
import styles from "./atlasModelChoice.module.css";
import { useHistory} from 'react-router-dom';
import { useEffect, useState } from 'react';
import ProjectMock from "shared/services/mock/projects"
import atlasPng from "assets/previewImages/atlas.png"
import ModelService from 'shared/services/Model.service';

function AtlasModelChoice({ 
    activeStep, setActiveStep, 
    selectedAtlas, setSelectedAtlas, 
    selectedModel, setSelectedModel, steps, path
}) {
    let [atlases, setAtlases] = useState([]);
    let [models, setModels] = useState([]);
    const history = useHistory();

    let headerStyle = {
        color: "#003560",
        fontSize: "1.6rem",
        fontWeight: "bold",
    }
    
    useEffect(() => {
        ProjectMock.getAtlases().then(a => setAtlases(a));
        ProjectMock.getModels().then(m => setModels(m));
    });

    return (
        <div>
            <Typography sx={headerStyle}>
                Pick an Atlas <HelpIcon sx={{color:"#B1CBDE"}} />
            </Typography>
            <hr className={styles.line}/>

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
                                imgLink={atlasPng}
                                selected={selectedAtlas===a.name}
                                onSelect={setSelectedAtlas}
                            />      
                        </Grid>
                    )
                }
            </Grid>
            
            <Typography sx={headerStyle} marginTop="20px">
                Pick a Model <HelpIcon sx={{color:"#B1CBDE"}} />
            </Typography>
            <hr className={styles.line}/>

            <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
            {
                models.map(m => 
                    <Grid item height="150px">
                    <ModelCardSelect width="225px"
                            height="97%"
                            title={m.name} 
                            description={m.description}
                            selected={selectedModel===m.name}
                            onSelect={setSelectedModel}
                    />               
                </Grid>)
            }
            </Grid>
            <Stack direction='row' justifyContent='space-between' sx={{ marginTop:'75px'}}>
                <CustomButton type='tertiary' onClick={() => history.push(`${path}`)}>
                    Cancel
                </CustomButton>
                <CustomButton type='primary' disabled={!selectedAtlas || !selectedModel} onClick={() => setActiveStep(1)}>
                    Confirm
                </CustomButton>
            </Stack>
        </div>
    )
}

export default AtlasModelChoice