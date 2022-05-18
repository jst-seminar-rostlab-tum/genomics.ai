import { Box, Stack } from "@mui/material";
import Typography from '@mui/material/Typography';
import { colors } from 'shared/theme/colors';

/** 
 * TabCard for TabGroup Component / List of cards with select feature used in FileUpload page
 * @param width
 * @param height
 * @data object containing information to be displayed on the card
 * @handleOnClick executed function on card click
 * @selected boolean parameter to indicate if a card element is selected
*/
export const TabCard = ({ width, height, data, handleOnClick, selected}) => {

    return (
        <Box onClick={handleOnClick} sx={{ width: {width}, height: {height}, backgroundColor: 'white', borderRadius: "0.625rem", marginTop: '0.5em', cursor: 'pointer', }}>
            <Box 
                sx={ selected ? {
                    boxShadow: "0px 0px 2px rgba(0,0,0, 0.15)", 
                    p: "0.5em", 
                    borderRadius: "0.625rem",
                    backgroundColor: colors.primary[300],
                    color: 'white'
                }
                : { 
                boxShadow: "0px 0px 2px rgba(0,0,0, 0.15)", 
                p: "0.5em", 
                borderRadius: "0.625rem",
                backgroundColor: 'white',
                ":hover": {
                    color: colors.primary,
                    ':hover': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
                    ':focus': { backgroundColor: colors.primary[300], transition: '0.4s', color: 'white' },
                    ':disabled': { backgroundColor: '#EBEFFF', transition: '0.4s', color: colors.primary[600] }
                }}}
            >
                <Stack direction="row" spacing={4} p='0.1em' pl='0.3em'>
                    <Typography variant='body1'>
                        {data.name}
                    </Typography>
                </Stack>
            </Box>
        </Box>
    );
}