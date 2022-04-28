import { Box, Button, Typography } from '@mui/material'
import DashboardHeader from 'components/DashboardHeader'
import Tag from 'components/Tag'

export default function TestPage(){
    return (
        <DashboardHeader title="Title" fontSize="4.2em" width="100%" height="200px">

            <Typography>Technische Universität München</Typography>
            
            <Box sx={{position: "relative", width: "100px", height: "30px", marginLeft: "10px"}}>
                <Tag content="Public" variant="primary-default" />
            </Box>

            <Button sx={{marginLeft: "73%"}}>Add</Button>
        </DashboardHeader>
    )
}