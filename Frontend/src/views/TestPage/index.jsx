import { Box } from '@mui/material'
import DashboardHeader from 'components/DashboardHeader'

export default function TestPage(){
    return (
        <Box sx={{width: "500px", height: "200px"}}>
            <DashboardHeader title="Title" fontSize="2.2em"></DashboardHeader>
        </Box>
    )
}