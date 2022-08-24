import React from 'react';
import { Box, Typography, Link } from '@mui/material';
import { colors } from 'shared/theme/colors';

export default function NotFound(){
    
    return (
        <Box sx={{position: "absolute", width: "100%", margin: 'auto', top: "30%", textAlign: "center"}}>
            <Typography fontSize="44px" fontWeight="bold">Error 404</Typography>
            <Typography fontSize="24px" fontWeight="bold">Not Found</Typography>
            <Typography marginTop="20px">The resource could not be found on the server!</Typography>
            <Link href="#" style={{ textDecorationColor: colors.primary[400], textDecorationThickness: "1px" }}>
                <Typography sx={{color: colors.primary[400], textDecoration: "none", cursor: 'pointer', ":hover": {color: colors.primary[500]} }}>
                Go Back Home
                </Typography>
            </Link>
        </Box>
    );
}
