import AtlasCard from "components/Cards/AtlasCard";
import {ModelCard} from "components/Cards/ModelCard";
import Stepper from "@mui/material/Stepper";
import Step from "@mui/material/Step";
import Box from "@mui/material/Box";
import StepLabel from "@mui/material/StepLabel";
import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";
import HelpIcon from '@mui/icons-material/Help';
import React, {useState} from "react";
import styles from "./atlasModelChoice.module.css";

function AtlasModelChoice(props) {
    let steps = ["Pick Atlas and Model", "Choose File and Project details"];
    let [activeStep, setActiveStep] = useState(0);
    let [completed, setCompleted] = useState([]);

    let headerStyle = {
        color: "#003560",
        fontSize: "1.6rem",
        fontWeight: "bold"
    }
    
    return (
        <div>
            <br/>
            <Box width="500px" margin="auto">
                <Stepper activeStep={activeStep}>
                    {steps.map((labelText, index) => {
                        return (
                            <Step index={index}>
                                <StepLabel>
                                    {labelText}
                                </StepLabel>
                            </Step>
                        )
                    })}
                </Stepper>
            </Box>
            <br/>
            <br/>
            <Typography marginLeft="170px" sx={headerStyle}>
                Pick an Atlas <HelpIcon sx={{color:"#B1CBDE"}} />
            </Typography>
            <hr className={styles.line}/>

            <Grid container spacing={2} width="84%" marginLeft="150px">
                <Grid item xs={2.4}>
                    <AtlasCard width="225px"
                        height="300px"
                        title="Human - PBDC" 
                    />                
                </Grid>
                <Grid item xs={2.4}>
                    <AtlasCard width="225px"
                        height="300px"
                        title="Human - PBDC" 
                    />                
                </Grid>
                <Grid item xs={2.4}>
                    <AtlasCard width="225px"
                        height="300px"
                        title="Human - PBDC" 
                    />                
                </Grid>
                <Grid item xs={2.4}>
                    <AtlasCard width="225px"
                        height="300px"
                        title="Human - PBDC" 
                    />                
                </Grid>
                <Grid item xs={2.4}>
                    <AtlasCard width="225px"
                        height="300px"
                        title="Human - PBDC" 
                    />                
                </Grid>
            </Grid>
            
            <Typography marginLeft="170px" sx={headerStyle} marginTop="20px">
                Pick a Model <HelpIcon sx={{color:"#B1CBDE"}} />
            </Typography>
            <hr className={styles.line}/>

            <Grid container spacing={2} width="84%" marginLeft="150px">
                <Grid item xs={2.4}>
                    <ModelCard width="225px"
                            height="150px"
                            title="Model 1" 
                            description="Lorem ipsum dolor sit amet, consetetur sadipscing"
                    />               
                </Grid>
                <Grid item xs={2.4}>
                    <ModelCard width="225px"
                            height="150px"
                            title="Model 2" 
                            description="Lorem ipsum dolor sit amet, consetetur sadipscing"
                    />               
                </Grid>
                <Grid item xs={2.4}>
                    <ModelCard width="225px"
                            height="150px"
                            title="Model 3" 
                            description="Lorem ipsum dolor sit amet, consetetur sadipscing"
                    />               
                </Grid>
            </Grid>
        </div>
    )
}

export default AtlasModelChoice